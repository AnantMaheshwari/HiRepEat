"""HiRepEat pipeline: tandem compression, LTR detection, clustering, consensus, iterative collapse.
"""

import os
import subprocess
import time
from pathlib import Path

from engine.constants import FASTAN_SCAN_FLANK, PRUNE_GENOMES
from seq.info import aln_info, fasta_size, scaffold_order, seq_size
from utilities import log
from utilities.final_library import is_ltr, write_final_library
from utilities.path_operations import gml_safe_label
from utilities.timing import timed

_REPO = Path(__file__).parent.parent
FASTAN_BIN       = str(_REPO / 'external' / 'FASTAN'   / 'FasTAN')
FASTLTR_BIN      = str(_REPO / 'external' / 'FASTAN'   / 'FastLTR')
TACO_BIN         = str(_REPO / 'external' / 'alntools' / 'taco')
EXTRACT_LTRS_BIN = str(_REPO / 'seq'      / 'extract_ltrs')


# ---------------------------------------------------------------------------
# Round 0: tandem compression, detection, extraction
# ---------------------------------------------------------------------------

def _tandem_dir(output_dir):
    """Directory holding the tandem-compression intermediate files.
    """
    return Path(output_dir) / 'tandem'


def _prune_genome_snapshot(path, keep):
    """Delete a superseded `.1seq` genome snapshot."""
    if path is None:
        return
    target = Path(path)
    if str(target) in keep or target.suffix != '.1seq' or not target.exists():
        return
    freed = target.stat().st_size
    target.unlink()
    log.info(f"  [prune] {target.name} ({freed / 1e6:.0f} MB)")


def run_fastan(genome, output_dir, max_iterations=1, skip_converged_contigs=True,
               prune_genomes=False):
    """Run FasTAN + taco compress iteratively until no new tandems are found.
    """
    stem = Path(genome).stem
    tandem_dir = _tandem_dir(output_dir)
    tandem_dir.mkdir(parents=True, exist_ok=True)

    current_genome = genome
    taco_layers = []   # (taco_file, seq_file) per iteration
    iter_stats  = []
    iteration = 0
    original_size = fasta_size(genome)
    prev_active_file = None   # <stem>_active_iter<N>.bed, consumed by iter N+1

    while iteration < max_iterations:
        iteration += 1
        iter_tag = f"_taniter{iteration}"
        tan_aln = str(tandem_dir / f"{stem}{iter_tag}.tan.1aln")
        prefix  = str(tandem_dir / f"{stem}{iter_tag}")

        log.info(f"  FasTAN iteration {iteration} on {Path(current_genome).name} ...")
        t_iter = time.perf_counter()
        cmd = [FASTAN_BIN, '-va', f'-o{tan_aln}']
        active_pct = 100.0
        skipped_bp = 0
        if skip_converged_contigs and prev_active_file is not None:
            cmd.append(f'-r{prev_active_file}')
            prior_info = iter_stats[-1].get('aln_info') if iter_stats else None
            if prior_info and prior_info['total_bp'] > 0:
                active_pct = 100.0 * prior_info['active_bp'] / prior_info['total_bp']
                skipped_bp = prior_info['total_bp'] - prior_info['active_bp']
        cmd.append(current_genome)
        with timed(f"FasTAN iter {iteration}"):
            log.run_tool(cmd, f"fastan_iter{iteration}")

        aln_summary = aln_info(tan_aln)
        n_alns = aln_summary['alignments']
        if n_alns == 0:
            log.step(f"  FasTAN iteration {iteration}: no new tandem repeats found. Converged.")
            cur_size = seq_size(current_genome)
            iter_stats.append({
                'iter': iteration, 'alignments': 0,
                'genome_bp': cur_size,
                'compress_pct': 100.0 * (1 - cur_size / original_size) if original_size else 0.0,
                'time_s': time.perf_counter() - t_iter,
                'active_pct': active_pct,
                'skipped_bp': skipped_bp,
                'aln_info': aln_summary,
            })
            break
        log.info(f"  FasTAN iteration {iteration}: {n_alns:,} tandem alignments found.")

        with timed(f"taco compress iter {iteration}"):
            log.run_tool(
                [TACO_BIN, 'compress', '-o', prefix, tan_aln, current_genome],
                f"taco_compress_iter{iteration}",
            )

        iter_taco   = f"{prefix}.1taco"
        iter_seq    = f"{prefix}.1seq"
        iter_genome = iter_seq
        taco_layers.append((iter_taco, iter_seq))

        if skip_converged_contigs:
            next_active = str(tandem_dir / f"{stem}_active_iter{iteration}.bed")
            with open(next_active, 'w') as bed:
                log.run_tool(
                    [TACO_BIN, 'windows', '-f', str(FASTAN_SCAN_FLANK), iter_taco],
                    f"taco_windows_iter{iteration}", stdout=bed,
                )
            prev_active_file = next_active

        prev_size = seq_size(current_genome)
        new_size  = seq_size(iter_genome)
        cumulative_pct = 100.0 * (1 - new_size / original_size) if original_size else 0.0
        iter_pct = 100.0 * (1 - new_size / prev_size) if prev_size else 0.0
        log.step(f"  FasTAN iteration {iteration}: {n_alns:,} tandems, "
                 f"{prev_size:,} → {new_size:,} bp ({iter_pct:.2f}% reduction)")

        iter_stats.append({
            'iter': iteration, 'alignments': n_alns,
            'genome_bp': new_size, 'compress_pct': cumulative_pct,
            'iter_pct': iter_pct,
            'time_s': time.perf_counter() - t_iter,
            'active_pct': active_pct,
            'skipped_bp': skipped_bp,
            'aln_info': aln_summary,
        })

        if prune_genomes and len(taco_layers) >= 2:
            _prune_genome_snapshot(taco_layers[-2][1], keep={genome})
            taco_layers[-2] = (taco_layers[-2][0], None)
        current_genome = iter_genome
    else:
        if max_iterations > 1:
            log.warn(f"FasTAN did not converge after {max_iterations} iterations.")

    log.info(f"\n{'=' * 72}")
    log.info("  Iterative FasTAN tandem compression summary")
    log.info(f"{'=' * 72}")
    log.info(f"  {'Iter':<6} {'Tandems':>8} {'Genome_Mbp':>11} "
             f"{'Iter%':>7} {'Cumul%':>7} {'Active%':>8} {'SkipMbp':>8} {'Time(s)':>8}")
    log.info(f"  {'-'*6} {'-'*8} {'-'*11} {'-'*7} {'-'*7} {'-'*8} {'-'*8} {'-'*8}")
    log.info(f"  {'orig':<6} {'-':>8} {original_size/1e6:>11.1f} "
             f"{'-':>7} {'-':>7} {'-':>8} {'-':>8} {'-':>8}")
    for s in iter_stats:
        alns_s = f"{s['alignments']:,}" if s['alignments'] > 0 else '0'
        gbp_s  = f"{s['genome_bp']/1e6:.1f}"
        ipct   = f"{s.get('iter_pct', 0.0):.1f}" if s['alignments'] > 0 else '-'
        cpct   = f"{s['compress_pct']:.1f}"
        apct   = f"{s.get('active_pct', 100.0):.1f}"
        skmbp  = f"{s.get('skipped_bp', 0)/1e6:.1f}"
        t_s    = s['time_s']
        log.info(f"  {s['iter']:<6} {alns_s:>8} {gbp_s:>11} "
                 f"{ipct:>7} {cpct:>7} {apct:>8} {skmbp:>8} {t_s:>8.1f}")
    log.info('')

    final_pct = iter_stats[-1]['compress_pct'] if iter_stats else 0.0
    log.stage(f"Tandem compression done: {iteration} iteration(s), "
              f"{len(taco_layers)} taco layer(s), {final_pct:.1f}% removed")
    log.info(f"  Tandem-compressed genome: {current_genome}")
    return current_genome, taco_layers


def run_fastltr(genome, output):
    """Run FastLTR on genome."""
    log.step(f"Running FastLTR on {Path(genome).name} ...")
    with timed("FastLTR"):
        log.run_tool([FASTLTR_BIN, '-v', genome, output], "fastltr_round1")


class ResumeError(RuntimeError):
    """A cached output directory cannot be safely resumed from."""


def check_resume(prior_params, genome, base, output_dir):
    """Some engineering to handle intermediately stopped runs, to pick up
    from already kept work (TODO: (anant): stress-test this in various error modes)
    """
    if genome is None:
        return                      
    if not (Path(output_dir) / 'round1' / f"{base}_ltr_sequences.fasta").exists():
        return
    _check_resume_genome(prior_params, genome, output_dir)
    _recover_tandem_layers(output_dir, Path(genome).stem)


def _check_resume_genome(prior_params, genome, output_dir):
    """Refuse to resume onto a different genome than the cache was built from.
    """
    prior_input = (prior_params or {}).get('input') or {}
    prior_genome = prior_input.get('genome')
    if not prior_genome:
        return                      # no record (or --genome-list): nothing to check

    now = str(Path(genome).resolve())
    prior_bp = prior_input.get('genome_bp')
    now_bp = seq_size(genome) if prior_bp is not None else None

    if now == prior_genome and (prior_bp is None or now_bp == prior_bp):
        return

    raise ResumeError(
        f"{output_dir} holds detections from a different genome:\n"
        f"    recorded: {prior_genome}"
        + (f" ({prior_bp:,} bp)" if isinstance(prior_bp, int) else '') + "\n"
        f"    now:      {now}"
        + (f" ({now_bp:,} bp)" if isinstance(now_bp, int) else '') + "\n"
        "Use a different output directory, or delete the existing one. Reusing "
        "it would lift this genome's coordinates through the other's taco stack."
    )


def _recover_tandem_layers(output_dir, stem):
    """Recover a previous run's taco layers, for resume.
    """
    layers = []
    i = 0
    while True:
        i += 1
        prefix = str(_tandem_dir(output_dir) / f"{stem}_taniter{i}")
        taco = f"{prefix}.1taco"
        if not Path(taco).exists():
            break
        seq = f"{prefix}.1seq"
        layers.append((taco, seq if Path(seq).exists() else None))

    if not layers:
        return [], None

    genome_for_ltr = layers[-1][1]
    if genome_for_ltr is None:
        raise ResumeError(
            f"{layers[-1][0]} is the last tandem layer but its "
            f"{Path(layers[-1][0]).with_suffix('.1seq').name} is missing — most "
            "likely pruned by a previous run with constants.PRUNE_GENOMES on.\n"
            "That file is the genome every later stage works in, and it cannot "
            "be reconstructed from what is left. Delete "
            f"{Path(output_dir).name} and re-run detection."
        )
    return layers, genome_for_ltr


def detect_and_extract(genome, base, output_dir, tandem_rounds,
                       skip_converged_contigs=True, prune_genomes=False,
                       detect_dir=None):
    """Run detection + extraction for one genome.
    """
    detect_dir = Path(output_dir if detect_dir is None else detect_dir)
    detect_dir.mkdir(parents=True, exist_ok=True)
    extracted_fasta = detect_dir / f"{base}_ltr_sequences.fasta"
    stem = Path(genome).stem

    if extracted_fasta.exists():
        log.stage(f"[resume] Reusing existing detections: {extracted_fasta}")
        log.step("         FasTAN, FastLTR and extraction are all skipped. "
                 "Delete that file to force re-detection.")
        tandem_layers, recovered = _recover_tandem_layers(output_dir, stem)
        if tandem_rounds > 0 and not tandem_layers:
            raise ResumeError(
                f"{extracted_fasta.name} exists, so a previous run reached "
                f"extraction, but no `{stem}_taniter*.1taco` layers were found "
                f"in {_tandem_dir(output_dir)}.\n"
            )
        genome_for_ltr = recovered if recovered is not None else genome
        log.step(f"         Recovered {len(tandem_layers)} taco layer(s); "
                 f"working genome {Path(genome_for_ltr).name}")
        return extracted_fasta, genome_for_ltr, tandem_layers

    t_total = time.perf_counter()

    tandem_layers = []
    genome_for_ltr = genome
    if tandem_rounds > 0:
        genome_for_ltr, tandem_layers = run_fastan(
            genome, output_dir,
            max_iterations=tandem_rounds,
            skip_converged_contigs=skip_converged_contigs,
            prune_genomes=prune_genomes,
        )

    ltr_file = detect_dir / f"{base}_LTR.1aln"
    run_fastltr(genome_for_ltr, str(ltr_file))

    with timed("extract TE sequences from reference"):
        log.run_tool(
            [EXTRACT_LTRS_BIN, genome_for_ltr, str(ltr_file), str(extracted_fasta)],
            "extract_ltrs_round1",
        )
    log.step(f"Extracted sequences to {extracted_fasta}")

    log.info(f"  [timing] detect_and_extract total: "
             f"{time.perf_counter() - t_total:.2f}s")
    return extracted_fasta, genome_for_ltr, tandem_layers


def guard_tandem_only_outputs(output_dir, genome, detect_marker):
    stem = Path(genome).stem
    if os.environ.get('HIREPEAT_FORCE_TANDEM'):
        log.warn("HIREPEAT_FORCE_TANDEM set: overwriting existing tandem layers.")
        return

    detect_marker = Path(detect_marker)
    if detect_marker.exists():
        raise ResumeError(
            f"{detect_marker} already exists"
        )

    existing = sorted(_tandem_dir(output_dir).glob(f"{stem}_taniter*.1taco"))
    if existing:
        raise ResumeError(
            f"{_tandem_dir(output_dir)} already holds {len(existing)} tandem layers"
        )


def run_tandem_only(genome, output_dir, tandem_rounds,
                    skip_converged_contigs=True, prune_genomes=False):
    """Tandem compression only — no LTR-RT detection.  Used by `--rounds 0`.
    """
    log.stage(f"Tandem compression only ({tandem_rounds} max iteration(s))")
    collapsed, layers = run_fastan(
        genome, output_dir,
        max_iterations=tandem_rounds,
        skip_converged_contigs=skip_converged_contigs,
        prune_genomes=prune_genomes,
    )
    if not layers:
        log.warn("FasTAN found no tandem repeats.")
    log.step(f"  Tandem-compressed genome: {collapsed}")
    log.step(f"  {len(layers)} taco layer(s) under {_tandem_dir(output_dir)}")
    return collapsed, layers


def combine_genome_fastas(extracted_fastas, combined_fasta):
    """Concatenate per-genome extracted FASTAs, tagging every record with its
    source genome (``<id>_<label>`` plus a ``genome=`` header field)."""
    with open(combined_fasta, 'w') as out:
        for fa, label in extracted_fastas:
            gml_label = gml_safe_label(label)
            with open(fa) as inp:
                for line in inp:
                    if line.startswith('>'):
                        seqid, _, rest = line[1:].partition(' ')
                        new_rest = (f"genome={gml_label} {rest}".strip()
                                    if rest.strip() else f"genome={gml_label}")
                        out.write(f'>{seqid}_{label} {new_rest}\n')
                    else:
                        out.write(line)


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_classifications(classified_fasta):
    """Parse consensus_name → classification from a classified FASTA.
    Handles the ``class=LTR/Gypsy`` key=value format produced by
    classify_tesorter.
    """
    classifications = {}
    with open(classified_fasta) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            parts = line[1:].strip().split()
            name = parts[0]
            cls = 'unknown'
            for p in parts[1:]:
                if p.startswith('class='):
                    cls = p.split('=', 1)[1]
                    break
                if p.startswith('#'):
                    cls = p.lstrip('#')
                    break
            classifications[name] = cls
    return classifications


def parse_cluster_membership(tsv_path):
    """Parse consensus_name → [member_ids] from cluster_membership.tsv."""
    membership = {}
    with open(tsv_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            membership[parts[0]] = parts[1].split(',')
    return membership


def parse_element_coordinates(extracted_fasta):
    """Parse element coordinates from structural FASTA headers.
    """
    coords = {}
    with open(extracted_fasta) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            parts = line[1:].strip().split()
            name = parts[0]
            meta = {}
            for p in parts[1:]:
                if '=' in p:
                    k, v = p.split('=', 1)
                    meta[k] = v

            scaffold = meta.get('scaffold', '')
            repeat1  = meta.get('repeat1', '')
            repeat2  = meta.get('repeat2', '')
            if not repeat1 or not repeat2:
                continue

            r1_start = int(repeat1.split('-')[0])
            r2_end   = int(repeat2.split('-')[1])
            coords[name] = (scaffold, r1_start, r2_end)
    return coords


def _merge_events(events):
    events.sort()
    merged = []
    for ev in events:
        if merged and ev[0] == merged[-1][0] and ev[1] <= merged[-1][2]:
            merged[-1] = (ev[0], merged[-1][1], max(ev[2], merged[-1][2]), 0)
        else:
            merged.append(ev)
    return merged


# ---------------------------------------------------------------------------
# Collapse
# ---------------------------------------------------------------------------

def write_detection_bed(extracted_fasta, genome_fasta, output_bed):
    """Write detected elements as BED for taco lift.
    """
    coords = parse_element_coordinates(extracted_fasta)
    scaffold_to_idx = {
        name: i for i, name in enumerate(scaffold_order(genome_fasta))
    }
    written = 0
    with open(output_bed, 'w') as f:
        for elem_id, (scaffold, start, end) in sorted(coords.items()):
            if scaffold not in scaffold_to_idx:
                continue
            f.write(f"{scaffold_to_idx[scaffold]}\t{start}\t{end}\t{elem_id}\n")
            written += 1
    return written


def write_flank_bed(extracted_fasta, genome_fasta, output_bed):
    """Write the flanks of each detection as BED for taco lift.
    """
    scaffold_to_idx = {
        name: i for i, name in enumerate(scaffold_order(genome_fasta))
    }
    written = 0
    with open(extracted_fasta) as f, open(output_bed, 'w') as out:
        for line in f:
            if not line.startswith('>'):
                continue
            parts = line[1:].strip().split()
            name = parts[0]
            meta = {}
            for p in parts[1:]:
                if '=' in p:
                    k, v = p.split('=', 1)
                    meta[k] = v
            scaffold = meta.get('scaffold', '')
            flank_l  = meta.get('flank_l', '')
            flank_r  = meta.get('flank_r', '')
            if scaffold not in scaffold_to_idx or not flank_l or not flank_r:
                continue
            foot_start = int(flank_l.split('-')[0])
            foot_end   = int(flank_r.split('-')[1])
            out.write(f"{scaffold_to_idx[scaffold]}\t{foot_start}\t{foot_end}\t{name}\n")
            written += 1
    return written


def _parse_lifted_bed(bed_path):
    """Parse a (possibly lifted) detection BED into element_id -> (seq_index, start, end).
    """
    out = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            out[parts[3]] = (int(parts[0]), int(parts[1]), int(parts[2]))
    return out


# ---------------------------------------------------------------------------
# Tool wrappers
# ---------------------------------------------------------------------------

def run_collapse(genome_fasta, spec_tsv, output_prefix, log_name='taco_collapse'):
    """Run taco collapse.  Returns (seq_file, taco_file).
    """
    cmd = [TACO_BIN, 'collapse', '-o', output_prefix, spec_tsv, genome_fasta]
    log.run_tool(cmd, log_name)
    return f"{output_prefix}.1seq", f"{output_prefix}.1taco"


def run_lift(annotations_bed, taco_stack, output_bed, log_name='taco_lift'):
    """Lift coordinates from collapsed space back to original via taco stack.
    """
    cmd = [TACO_BIN, 'lift', annotations_bed] + list(taco_stack)
    with open(output_bed, 'w') as out:
        log.run_tool(cmd, log_name, stdout=out)


def run_info(taco_file):
    """Run taco info and return output as string."""
    result = subprocess.run(
        [TACO_BIN, 'info', taco_file],
        capture_output=True, text=True, check=True,
    )
    return result.stdout.strip()


# ---------------------------------------------------------------------------
# Cross-round accumulation helpers
# ---------------------------------------------------------------------------

_is_ltr_class = is_ltr


def _read_fasta_records(path):
    """Read a FASTA into an ordered dict: name -> (header_no_newline, seq)."""
    out = {}
    name = header = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out[name] = (header, ''.join(buf))
                header = line.rstrip('\n')
                name   = header[1:].split()[0]
                buf    = []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = (header, ''.join(buf))
    return out


def _collapse_spec_from_coords(member_coords, genome_fasta, output_tsv):
    """Write a taco collapse spec from explicit (scaffold, start, end) tuples.
    """
    order           = scaffold_order(genome_fasta)
    scaffold_to_idx = {name: i for i, name in enumerate(order)}

    events = []
    for scaffold, start, end in member_coords:
        if scaffold in scaffold_to_idx:
            events.append((scaffold_to_idx[scaffold], start, end, 0))

    merged = _merge_events(events)
    with open(output_tsv, 'w') as f:
        f.write("# seq_index\torig_start\torig_end\tunit_len\n")
        for seq_idx, start, end, unit in merged:
            f.write(f"{seq_idx}\t{start}\t{end}\t{unit}\n")

    active = {order[ev[0]] for ev in merged if 0 <= ev[0] < len(order)}
    return len(merged), active


# ---------------------------------------------------------------------------
# Main iterative loop
# ---------------------------------------------------------------------------

def build_feature_config(args):
    """FeatureClusterConfig for the clustering stage.
    """
    from engine.feature_cluster import FeatureClusterConfig
    return FeatureClusterConfig(
        threads=args.threads,
        keep_alignments=args.keep_alignments,
        debug=args.debug,
        use_flank_filter=args.flank_filter,
    )


def classify(consensi_fasta, output_fasta, args, label='classify',
             log_name='tesorter'):
    """Classify a consensus library with TEsorter.  Returns the output path.
    """
    from engine.classify_tesorter import classify_consensi_tesorter
    from engine.constants import TESORTER_DB, tesorter_gydb
    # The primary database first; the gydb pass appends a second run whose calls
    # are only ever consulted where the primary found no LTR family.
    dbs = [TESORTER_DB]
    if tesorter_gydb():
        dbs.append('gydb')
    with timed(f"{label} (tesorter)"):
        classify_consensi_tesorter(
            input_fasta=consensi_fasta, output_fasta=output_fasta,
            db=dbs, threads=args.threads, log_name=log_name,
        )
    return output_fasta


def cluster_and_classify(input_fasta, output_dir, args):
    """Cluster one FASTA into families and classify the consensi.
    """
    os.environ['OMP_NUM_THREADS'] = str(max(1, args.threads))

    from engine.feature_cluster import (
        _load_records, cluster_records_incremental, _ALGO_DIR, _write_gml,
        write_classified_gml,
    )
    from algorithms.c_libs import KmerLib, WFALib

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fcfg = build_feature_config(args)

    with timed("load FASTA + parse structural coordinates"):
        records = _load_records(input_fasta)
    if not records:
        empty = str(out / 'consensi.fasta')
        open(empty, 'w').close()
        return empty, empty

    kmer_lib = KmerLib(_ALGO_DIR / 'kmer_seed.so')
    wfa_lib  = WFALib(_ALGO_DIR / 'wfa_align.so')

    with timed("cluster"):
        consensi_fa, membership_tsv, _all_results, G = cluster_records_incremental(
            records, fcfg, kmer_lib, wfa_lib, out,
        )

    if fcfg.debug:
        _write_gml(G, out)

    classified_fa = str(out / 'consensi_classified_tesorter.fasta')
    if Path(consensi_fa).exists() and Path(consensi_fa).stat().st_size > 0:
        classify(consensi_fa, classified_fa, args)
        if fcfg.debug and (out / 'similarity_graph.gml').exists():
            write_classified_gml(
                gml_path=str(out / 'similarity_graph.gml'),
                membership_tsv=membership_tsv,
                classified_fasta=classified_fa,
                output_path=str(out / 'similarity_graph_classified.gml'),
            )
    else:
        open(classified_fa, 'w').close()

    return consensi_fa, classified_fa




def run_rounds(genome_for_ltr, base, output_dir, round1_extracted,
               original_genome_size, args, tandem_layers=None):
    """Exact incremental TE collapse via a persistent element similarity graph.
    """
    # libgomp reads OMP_NUM_THREADS exactly once, in its load-time constructor,
    # so this MUST run before anything for parallelism
    os.environ['OMP_NUM_THREADS'] = str(max(1, args.threads))

    from engine.feature_cluster import (
        _load_records, cluster_records_incremental, _ALGO_DIR,
        _write_gml, write_classified_gml,
    )
    from algorithms.c_libs import KmerLib, WFALib

    kmer_lib = KmerLib(_ALGO_DIR / 'kmer_seed.so')
    wfa_lib  = WFALib(_ALGO_DIR / 'wfa_align.so')
    fcfg     = build_feature_config(args)

    # ── persistent graph state across rounds ────────────────────────────────
    node_seq:        dict = {}    # name -> (header, seq) for every current node
    contig_nodes:    dict = {}    # scaffold -> set[name] currently in the genome
    collapsed_nodes: set  = set() # names collapsed out of the genome (frozen)
    reuse_store:     dict = {}    # (name_a, name_b) -> edge attrs (surviving edges)
    frozen_params:   dict = {}    # set-dependent seeding params, frozen at round 1
    node_orig_coords: dict = {}   # name -> (seq_index, orig_start, orig_end) in original space
    last_removed_orig: dict = {}  # seq_index -> [(ra, rb), ...] prev round's removals (original space)
    scaffold_count            = None 

    taco_stack     = []
    lifted_beds    = []
    round_stats    = []
    current_genome = str(genome_for_ltr)
   
    keep_genomes   = {str(genome_for_ltr), str(getattr(args, 'genome', '') or '')}

    skip_converged   = True
    next_active_file = None

    no_node_reuse = bool(os.environ.get('HIREPEAT_NO_NODE_REUSE'))
    if no_node_reuse:
        log.stage("  [HIREPEAT_NO_NODE_REUSE] node + edge reuse DISABLED "
                  "(full universe rebuild + from-scratch re-cluster each round)")

    final_classified_src = None
    final_membership_src = None

    # ── Round 0: tandem compression layers ──────────────────────────────────
    if tandem_layers:
        for taco_file, _ in tandem_layers:
            taco_stack.append(taco_file)
        compressed_size = seq_size(current_genome)
        compress_pct = 100.0 * (1 - compressed_size / original_genome_size)
        log.info(f"\n  Round 0 (tandem compression, {len(tandem_layers)} layer(s)):")
        try:
            log.block(log.INFO, run_info(tandem_layers[-1][0]))
        except Exception:
            pass
        log.info(f"  Genome size: {original_genome_size:,} → {compressed_size:,} bp "
                 f"({compress_pct:.1f}% removed)")
        round_stats.append({
            'round': '0(tan)', 'detections': '-', 'families': '-',
            'non_ltr': '-', 'new_ltr_members': '-', 'collapse_events': '-',
            'genome_size_bp': compressed_size, 'compression_pct': compress_pct,
            'time_s': 0.0,
        })

    max_rounds = args.rounds

    for round_num in range(1, max_rounds + 1):
        t_round = time.perf_counter()
        round_dir = output_dir / f"round{round_num}"
        round_dir.mkdir(exist_ok=True)

        stack_note = f", {len(taco_stack)} taco layer(s)" if taco_stack else ''
        log.stage(f"\n── Round {round_num}/{max_rounds} ── "
                  f"{Path(current_genome).name}{stack_note}")
        log.info(f"  Genome: {current_genome}")

        # ── detect + extract for this round ─────────────────────────────────
        round_active_pct = 100.0
        round_skipped_bp = 0
        if round_num == 1:
            extracted_fa = str(round1_extracted)
        else:
            ltr_file = str(round_dir / f"{base}_round{round_num}_LTR.1aln")
            if Path(ltr_file).exists():
                try:
                    n_prev = aln_info(ltr_file)['alignments']
                except Exception:
                    n_prev = None
                log.stage(f"  [resume] Reusing existing {Path(ltr_file).name}"
                          + (f" ({n_prev:,} alignments)" if n_prev is not None
                             else " (unreadable — delete it if this run fails)"))
                if n_prev == 0:
                    log.warn(f"{ltr_file} contains no alignments; delete it to "
                             f"re-run FastLTR for this round.")
            else:
                cmd = [FASTLTR_BIN, '-v']
                if skip_converged and next_active_file is not None \
                        and Path(next_active_file).exists():
                    cmd.append(f'-r{next_active_file}')
                    try:
                        with open(next_active_file) as af:
                            active_names = {l.strip() for l in af if l.strip()}
                        all_scafs = scaffold_order(current_genome)
                        if all_scafs:
                            round_active_pct = 100.0 * len(active_names) / len(all_scafs)
                        cur_size = seq_size(current_genome)
                        round_skipped_bp = int(cur_size * (1 - round_active_pct / 100.0))
                    except Exception:
                        pass
                cmd.extend([current_genome, ltr_file])
                log.step(f"  Running FastLTR on {Path(current_genome).name} ...")
                with timed(f"round {round_num}: FastLTR"):
                    log.run_tool(cmd, f"fastltr_round{round_num}")

            extracted_fa = str(round_dir / f"{base}_round{round_num}_ltr_sequences.fasta")
            with timed(f"round {round_num}: extract sequences"):
                log.run_tool(
                    [EXTRACT_LTRS_BIN, current_genome, ltr_file, extracted_fa],
                    f"extract_ltrs_round{round_num}",
                )

        n_new = sum(1 for line in open(extracted_fa) if line.startswith('>'))
        log.step(f"  {n_new:,} detections extracted in round {round_num}")

        # ── BED of this round's detections, lifted to original coords ───────
        round_bed = str(round_dir / f"{base}_round{round_num}_detections.bed")
        write_detection_bed(extracted_fa, current_genome, round_bed)
        if taco_stack:
            lifted_bed = str(round_dir / f"{base}_round{round_num}_detections_original.bed")
            run_lift(round_bed, taco_stack, lifted_bed,
                     log_name=f"taco_lift_round{round_num}_detections")
            lifted_beds.append(lifted_bed)
        else:
            lifted_beds.append(round_bed)

        orig_elem = _parse_lifted_bed(lifted_beds[-1]) 

        orig_flank: dict = {}
        if taco_stack:
            flank_bed = str(round_dir / f"{base}_round{round_num}_flank.bed")
            write_flank_bed(extracted_fa, current_genome, flank_bed)
            flank_lifted = str(round_dir / f"{base}_round{round_num}_flank_original.bed")
            run_lift(flank_bed, taco_stack, flank_lifted,
                     log_name=f"taco_lift_round{round_num}_flank")
            orig_flank = _parse_lifted_bed(flank_lifted)

        # ── reconcile graph nodes ───────────────────────────────────────────
        raw    = _read_fasta_records(extracted_fa)        # orig_name -> (header, seq)
        coords = parse_element_coordinates(extracted_fa)  # orig_name -> (scaf,s,e)
        coords_current: dict = {}    # node_name -> (scaf,s,e) for THIS round's detections
        prefix    = f"r{round_num}_"
        by_contig: dict = {}
        for orig_name, (header, seq) in raw.items():
            scaffold = coords.get(orig_name, ('', 0, 0))[0]
            if not scaffold:
                continue
            pname   = prefix + orig_name
            pheader = header.replace(f">{orig_name}", f">{pname}", 1)
            by_contig.setdefault(scaffold, {})[pname] = (pheader, seq)


        cur_order = scaffold_order(current_genome)
        if scaffold_count is None:
            scaffold_count = len(cur_order)
        elif len(cur_order) != scaffold_count:
            raise RuntimeError(
                f"scaffold count changed {scaffold_count} → {len(cur_order)} at "
                f"round {round_num}; seq_index reuse keys would be invalid")


        if skip_converged and next_active_file and Path(next_active_file).exists():
            with open(next_active_file) as af:
                reran_contigs = {l.strip() for l in af if l.strip()}
        else:
            reran_contigs = set(cur_order)

        new_names: set = set()
        if round_num == 1 or not taco_stack or no_node_reuse:
            for scaffold in sorted(reran_contigs):
                for old in contig_nodes.get(scaffold, set()):
                    node_seq.pop(old, None)
                    node_orig_coords.pop(old, None)
                contig_nodes[scaffold] = set()
                for pname, (pheader, seq) in by_contig.get(scaffold, {}).items():
                    node_seq[pname] = (pheader, seq)
                    contig_nodes[scaffold].add(pname)
                    coords_current[pname] = coords[pname[len(prefix):]]
                    oc = orig_elem.get(pname[len(prefix):])
                    if oc is not None:
                        node_orig_coords[pname] = oc
                    new_names.add(pname)
        else:
            # Incremental reuse path
            n_reused = 0
            for scaffold in sorted(reran_contigs):
                old_by_coord: dict = {}
                for old in sorted(contig_nodes.get(scaffold, set())):
                    oc = node_orig_coords.get(old)
                    if oc is not None:
                        old_by_coord.setdefault(oc, []).append(old)
                survivors: set = set()
                for pname, (pheader, seq) in by_contig.get(scaffold, {}).items():
                    orig_name = pname[len(prefix):]
                    es = orig_elem.get(orig_name)   
                    fs = orig_flank.get(orig_name)   
                    affected = False
                    if fs is not None:
                        si, fa, fb = fs
                        for ra, rb in last_removed_orig.get(si, ()):  # half-open
                            if fa < rb and ra < fb:
                                affected = True
                                break
                    match = None
                    if not affected and es is not None:
                        cands = old_by_coord.get(es)
                        if cands:
                            for i, cand in enumerate(cands):
                                if node_seq.get(cand, (None, None))[1] == seq:
                                    match = cand
                                    cands.pop(i)
                                    break
                    if match is not None:
                        coords_current[match] = coords[orig_name]
                        survivors.add(match)
                        n_reused += 1
                    else:
                        node_seq[pname] = (pheader, seq)
                        coords_current[pname] = coords[orig_name]
                        if es is not None:
                            node_orig_coords[pname] = es
                        new_names.add(pname)
                        survivors.add(pname)
                # drop old nodes neither reused nor re-detected on this contig
                for old in contig_nodes.get(scaffold, set()):
                    if old not in survivors:
                        node_seq.pop(old, None)
                        node_orig_coords.pop(old, None)
                contig_nodes[scaffold] = survivors
            log.info(f"  [incremental] reused {n_reused} node(s), {len(new_names)} new "
                     f"across {len(reran_contigs)} active contig(s)")

        # ── build the current universe of elements (stable + collapsed-back + new) ──────
        universe_fa = str(round_dir / f"{base}_round{round_num}_universe.fasta")
        with open(universe_fa, 'w') as out:
            for name, (header, seq) in node_seq.items():
                out.write(header + '\n')
                out.write(seq + '\n')
        records = _load_records(universe_fa)

        Path(universe_fa).unlink(missing_ok=True)
        n_stable = len(node_seq) - len(new_names) - len(collapsed_nodes)
        log.info(f"  Universe: {len(records)} nodes "
                 f"({len(new_names)} new, {len(collapsed_nodes)} collapsed-back, "
                 f"{n_stable} stable)")

        # ── incremental cluster (WFA only on new-incident pairs) ────────────
        reuse_edges = {(a, b): attrs for (a, b), attrs in reuse_store.items()
                       if a in node_seq and b in node_seq}
        first_round = (round_num == 1 or not reuse_edges or no_node_reuse)
        clustering_dir = round_dir / 'clustering'
        with timed(f"round {round_num}: incremental cluster"):
            consensi_fa, membership_tsv, all_results, G = cluster_records_incremental(
                records, fcfg, kmer_lib, wfa_lib, clustering_dir,
                new_names=(None if first_round else new_names),
                reuse_edges=(None if first_round else reuse_edges),
                frozen_params=frozen_params,
            )

        if fcfg.debug:
            _write_gml(G, clustering_dir)

        idx_name = [r.name for r in records]
        reuse_store = {
            (idx_name[u], idx_name[v]): attrs
            for u, v, attrs in G.edges(data=True)
        }

        # ── classify component consensi to decide what to collapse ──────────
        classifications = {}
        classified_fa = str(clustering_dir / 'consensi_classified_tesorter.fasta')
        if Path(consensi_fa).exists() and Path(consensi_fa).stat().st_size > 0:
            classify(consensi_fa, classified_fa, args,
                     label=f"round {round_num}: classify",
                     log_name=f"tesorter_round{round_num}")
            classifications = parse_classifications(classified_fa)
            # Debug-only: annotate the similarity graph with classifier labels.
            if fcfg.debug:
                gml_in = clustering_dir / 'similarity_graph.gml'
                if gml_in.exists():
                    write_classified_gml(
                        gml_path=str(gml_in),
                        membership_tsv=membership_tsv,
                        classified_fasta=classified_fa,
                        output_path=str(clustering_dir / 'similarity_graph_classified.gml'),
                    )
        membership = (parse_cluster_membership(membership_tsv)
                      if Path(membership_tsv).exists() else {})

        final_classified_src = classified_fa if Path(classified_fa).exists() else consensi_fa
        final_membership_src = membership_tsv

        n_ltr_fams = sum(1 for cls in classifications.values() if _is_ltr_class(cls))
        n_non_ltr  = len(classifications) - n_ltr_fams
        log.step(f"  Round {round_num}: {len(membership)} families "
                 f"({n_ltr_fams} LTR, {n_non_ltr} non-LTR)")

        # ── collapse this round's LTR-classified detections ─────────────────
        added_coords    = []
        collapsed_names: set = set()
        for cons_name, members in membership.items():
            if not _is_ltr_class(classifications.get(cons_name, 'unknown')):
                continue
            for m in members:
                if m in coords_current:        # only THIS round's fresh detections
                    added_coords.append(coords_current[m])
                    collapsed_names.add(m)

        n_events = 0
        last_removed_orig = {}  
        if round_num < max_rounds and added_coords:
            spec_tsv = str(round_dir / f"{base}_round{round_num}_collapse.tsv")
            output_prefix = str(round_dir / f"{base}_round{round_num}")
            with timed(f"round {round_num}: build collapse spec"):
                n_events, active_scaffolds = _collapse_spec_from_coords(
                    added_coords, current_genome, spec_tsv,
                )
            next_active_file = str(round_dir / f"active_round{round_num + 1}.txt")
            with open(next_active_file, 'w') as af:
                for name in sorted(active_scaffolds):
                    af.write(name + '\n')

            if n_events > 0:
                # Lift this round's removed intervals to original space for the next round.
                if taco_stack:
                    removed_lifted = str(
                        round_dir / f"{base}_round{round_num}_removed_original.bed")
                    run_lift(spec_tsv, taco_stack, removed_lifted,
                             log_name=f"taco_lift_round{round_num}_removed")
                    removed_src = removed_lifted
                else:
                    removed_src = spec_tsv   # current space == original space
                with open(removed_src) as rf:
                    for line in rf:
                        if line.startswith('#'):
                            continue
                        parts = line.split('\t')
                        if len(parts) < 3:
                            continue
                        si, ra, rb = int(parts[0]), int(parts[1]), int(parts[2])
                        last_removed_orig.setdefault(si, []).append((ra, rb))

                with timed(f"round {round_num}: taco collapse"):
                    collapsed_fa, taco_file = run_collapse(
                        current_genome, spec_tsv, output_prefix,
                        log_name=f"taco_collapse_round{round_num}",
                    )
                taco_stack.append(taco_file)
                if PRUNE_GENOMES:
                    _prune_genome_snapshot(current_genome, keep=keep_genomes)
                current_genome = collapsed_fa
                new_size = seq_size(collapsed_fa)
                log.step(f"  Collapsed {n_events:,} interval(s) → genome {new_size:,} bp "
                         f"({100.0 * (1 - new_size / original_genome_size):.1f}% removed)")

                # Collapsed detections become frozen nodes: out of the genome,
                # but kept as clustering evidence across future rounds.
                for m in collapsed_names:
                    scaf = coords_current[m][0]
                    contig_nodes.get(scaf, set()).discard(m)
                    collapsed_nodes.add(m)

        genome_size_now = seq_size(current_genome)
        round_stats.append({
            'round': round_num,
            'detections': n_new,
            'families': n_ltr_fams,
            'non_ltr': n_non_ltr,
            'new_ltr_members': len(collapsed_names),
            'collapse_events': n_events,
            'matched': '-',
            'genome_size_bp': genome_size_now,
            'compression_pct': 100.0 * (1 - genome_size_now / original_genome_size),
            'active_pct': round_active_pct,
            'skipped_bp': round_skipped_bp,
            'time_s': time.perf_counter() - t_round,
        })

        if round_num < max_rounds and not added_coords:
            log.stage("  No new LTR members to collapse — converged.")
            break

    # ── combine per-round lifted BEDs ───────────────────────────────────────
    all_rounds_dir = Path(output_dir) / 'all_rounds'
    all_rounds_dir.mkdir(parents=True, exist_ok=True)

    combined_bed = str(all_rounds_dir / f"{base}_all_rounds_original.bed")
    with open(combined_bed, 'w') as out:
        out.write("#seq_index\tstart\tend\telement_id\tround\n")
        for rnd, bed in enumerate(lifted_beds, start=1):
            if Path(bed).exists():
                with open(bed) as inp:
                    for line in inp:
                        if line.startswith('#'):
                            continue
                        out.write(f"{line.rstrip()}\t{rnd}\n")

    # ── final library = last round's consensi (deduplicated by construction) ─
    final_classified = str(all_rounds_dir / f"{base}_all_rounds_classified.fasta")
    final_membership = str(all_rounds_dir / f"{base}_all_rounds_membership.tsv")
    if final_classified_src and Path(final_classified_src).exists():
        Path(final_classified).write_bytes(Path(final_classified_src).read_bytes())
    else:
        open(final_classified, 'w').close()
    if final_membership_src and Path(final_membership_src).exists():
        Path(final_membership).write_bytes(Path(final_membership_src).read_bytes())
    else:
        open(final_membership, 'w').close()
    log.info(f"\n  all_rounds/ consensi:     {final_classified}")
    log.info(f"  all_rounds/ membership:   {final_membership}")
    log.info(f"  all_rounds/ detections:   {combined_bed}")

    ltr_fa, rm_fa, n_kept, n_dropped = write_final_library(
        final_classified, output_dir, base)
    log.info(f"  {n_dropped} consensi dropped as non-LTR, {n_kept} kept")
    log.stage(f"\n  {ltr_fa} is the final output LTR-RT library "
              f"({n_kept} families)")
    log.stage(f"  {rm_fa}   (RepeatMasker: >NAME#CLASS)") # writing an additional file compatible with RepeatMasker

    _print_incremental_summary(round_stats)

    return [final_classified], [final_membership], taco_stack, round_stats


def _print_incremental_summary(round_stats):
    """Per-round summary table for the incremental loop."""
    log.stage(f"\n{'=' * 80}")
    log.stage("  Incremental collapse summary")
    log.stage(f"{'=' * 80}")
    log.stage(f"  {'Round':<9} {'Detect':>7} {'Matched':>8} {'LTRfam':>7} "
              f"{'non-LTR':>8} {'NewMem':>7} {'Collapsed':>10} {'Genome_Mbp':>11} "
              f"{'Active%':>8} {'Time(s)':>8}")
    log.stage(f"  {'-'*9} {'-'*7} {'-'*8} {'-'*7} {'-'*8} {'-'*7} {'-'*10} "
              f"{'-'*11} {'-'*8} {'-'*8}")
    for s in round_stats:
        det   = s.get('detections', '-')
        mat   = s.get('matched', '-')
        ltr_  = s.get('families', '-')
        nlt   = s.get('non_ltr', '-')
        newm  = s.get('new_ltr_members', '-')
        coll  = s.get('collapse_events', '-')
        gbp   = s.get('genome_size_bp', 0)
        apct  = s.get('active_pct', None)
        t_s   = s.get('time_s', 0.0)
        gbp_s  = f"{gbp/1e6:.1f}" if isinstance(gbp, (int, float)) else '-'
        apct_s = f"{apct:.1f}" if isinstance(apct, float) else '-'
        log.stage(f"  {str(s['round']):<9} {str(det):>7} {str(mat):>8} {str(ltr_):>7} "
                  f"{str(nlt):>8} {str(newm):>7} {str(coll):>10} {gbp_s:>11} "
                  f"{apct_s:>8} {t_s:>8.1f}")
    log.stage('')
