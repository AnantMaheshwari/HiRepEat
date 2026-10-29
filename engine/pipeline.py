"""HiRepEat pipeline: tandem compression, LTR detection, clustering, consensus, iterative collapse.
"""

import os
import subprocess
import time
from dataclasses import dataclass, field
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


def check_genome_list_resume(prior_params, genome_paths, output_dir):
    """Refuse to resume a cross-genome run onto a different set of genomes.
    """
    prior = ((prior_params or {}).get('input') or {}).get('genomes')
    if not prior:
        return
    prior_set = {str(p) for p in (prior if isinstance(prior, list) else [prior])}
    now_set   = {str(Path(p).resolve()) for p in genome_paths}
    if prior_set == now_set:
        return

    added   = sorted(now_set - prior_set)
    removed = sorted(prior_set - now_set)
    raise ResumeError(
        f"{output_dir} holds results for a different set of genomes:\n"
        + ''.join(f"    added:   {p}\n" for p in added)
        + ''.join(f"    removed: {p}\n" for p in removed)
    )

def guard_tandem_only_outputs(output_dir, genome, detect_markers):
    stem = Path(genome).stem
    if os.environ.get('HIREPEAT_FORCE_TANDEM'):
        log.warn("HIREPEAT_FORCE_TANDEM set: overwriting existing tandem layers.")
        return

    if isinstance(detect_markers, (str, Path)):
        detect_markers = [detect_markers]
    for detect_marker in map(Path, detect_markers):
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
    n_dropped = 0
    for scaffold, start, end in member_coords:
        if scaffold in scaffold_to_idx:
            events.append((scaffold_to_idx[scaffold], start, end, 0))
        else:
            n_dropped += 1

    merged = _merge_events(events)
    with open(output_tsv, 'w') as f:
        f.write("# seq_index\torig_start\torig_end\tunit_len\n")
        for seq_idx, start, end, unit in merged:
            f.write(f"{seq_idx}\t{start}\t{end}\t{unit}\n")

    active = {order[ev[0]] for ev in merged if 0 <= ev[0] < len(order)}
    return len(merged), active, n_dropped


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


# ---------------------------------------------------------------------------
# Per-genome state for the iterative collapse loop
# ---------------------------------------------------------------------------

@dataclass
class GenomeState:
    """Collapse loop state across rounds.
    One instance in single-genome mode; one per genome under ``--genome-list``.
    """
    source_path:  str
    base:         str           
    original_size: int
    current_genome: str
    round1_extracted: str

    tandem_layers:     list = field(default_factory=list)
    taco_stack:        list = field(default_factory=list)
    node_seq:          dict = field(default_factory=dict)
    contig_nodes:      dict = field(default_factory=dict)
    collapsed_nodes:   set  = field(default_factory=set)
    node_orig_coords:  dict = field(default_factory=dict)
    last_removed_orig: dict = field(default_factory=dict)
    coords_current:    dict = field(default_factory=dict)
    scaffold_count:    int  = None
    next_active_file:  str  = None
    lifted_beds:       list = field(default_factory=list)
    converged:         bool = False

    # per-round scratch, consumed by the stats row
    n_new:            int   = 0
    round_active_pct: float = 100.0
    round_skipped_bp: int   = 0

    def round_dir(self, output_dir, round_num, cross_genome):
        d = Path(output_dir) / f"round{round_num}"
        return d / self.base if cross_genome else d

    def tool_tag(self, cross_genome):
        return f"_{self.base}" if cross_genome else ''

    def node_name(self, orig_name, round_num, cross_genome):
        """Graph-node name for a detection.
        """
        return (f"r{round_num}_{orig_name}_{self.base}" if cross_genome
                else f"r{round_num}_{orig_name}")

    def node_header(self, header, orig_name, pname, cross_genome):
        """Rewrite a detection's FASTA header to carry its node name.
        """
        out = header.replace(f">{orig_name}", f">{pname}", 1)
        if not cross_genome:
            return out
        seqid, _, rest = out[1:].partition(' ')
        tag = f"genome={gml_safe_label(self.base)}"
        return (f">{seqid} {tag} {rest.strip()}" if rest.strip()
                else f">{seqid} {tag}")


# ---------------------------------------------------------------------------
# One round, one genome: detect, lift, reconcile nodes
# ---------------------------------------------------------------------------

def _round0_tandem(st, tandem_layers, round_stats, cross_genome):
    """Seed a genome's taco stack with its tandem-compression layers."""
    if not tandem_layers:
        return
    for taco_file, _ in tandem_layers:
        st.taco_stack.append(taco_file)
    compressed_size = seq_size(st.current_genome)
    compress_pct = 100.0 * (1 - compressed_size / st.original_size)
    tag = f" [{st.base}]" if cross_genome else ''
    log.info(f"\n  Round 0 (tandem compression, {len(tandem_layers)} layer(s)):{tag}")
    try:
        log.block(log.INFO, run_info(tandem_layers[-1][0]))
    except Exception:
        pass
    log.info(f"  Genome size: {st.original_size:,} → {compressed_size:,} bp "
             f"({compress_pct:.1f}% removed)")
    row = {
        'round': '0(tan)', 'detections': '-', 'families': '-',
        'non_ltr': '-', 'new_ltr_members': '-', 'collapse_events': '-',
        'genome_size_bp': compressed_size, 'compression_pct': compress_pct,
        'time_s': 0.0,
    }
    if cross_genome:
        row['genome'] = st.base
    round_stats.append(row)


def _detect_round(st, round_num, max_rounds, round_dir, args,
                  skip_converged, cross_genome):
    st.round_active_pct = 100.0
    st.round_skipped_bp = 0

    stack_note = f", {len(st.taco_stack)} taco layer(s)" if st.taco_stack else ''
    log.stage(f"\n── Round {round_num}/{max_rounds} ── "
              f"{Path(st.current_genome).name}{stack_note}")
    log.info(f"  Genome: {st.current_genome}")

    if round_num == 1:
        extracted_fa = str(st.round1_extracted)
    else:
        round_dir.mkdir(parents=True, exist_ok=True)
        tag = st.tool_tag(cross_genome)
        ltr_file = str(round_dir / f"{st.base}_round{round_num}_LTR.1aln")
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
            if skip_converged and st.next_active_file is not None \
                    and Path(st.next_active_file).exists():
                cmd.append(f'-r{st.next_active_file}')
                try:
                    with open(st.next_active_file) as af:
                        active_names = {l.strip() for l in af if l.strip()}
                    all_scafs = scaffold_order(st.current_genome)
                    if all_scafs:
                        st.round_active_pct = 100.0 * len(active_names) / len(all_scafs)
                    cur_size = seq_size(st.current_genome)
                    st.round_skipped_bp = int(cur_size * (1 - st.round_active_pct / 100.0))
                except Exception:
                    pass
            cmd.extend([st.current_genome, ltr_file])
            log.step(f"  Running FastLTR on {Path(st.current_genome).name} ...")
            with timed(f"round {round_num}: FastLTR"):
                log.run_tool(cmd, f"fastltr_round{round_num}{tag}")

        extracted_fa = str(round_dir / f"{st.base}_round{round_num}_ltr_sequences.fasta")
        with timed(f"round {round_num}: extract sequences"):
            log.run_tool(
                [EXTRACT_LTRS_BIN, st.current_genome, ltr_file, extracted_fa],
                f"extract_ltrs_round{round_num}{tag}",
            )

    st.n_new = sum(1 for line in open(extracted_fa) if line.startswith('>'))
    log.step(f"  {st.n_new:,} detections extracted in round {round_num}")
    return extracted_fa


def _lift_round_beds(st, extracted_fa, round_num, round_dir, cross_genome):
    """Write this round's detection + flank and taco lift them to original space.
    """
    round_dir.mkdir(parents=True, exist_ok=True)
    tag = st.tool_tag(cross_genome)
    round_bed = str(round_dir / f"{st.base}_round{round_num}_detections.bed")
    write_detection_bed(extracted_fa, st.current_genome, round_bed)
    if st.taco_stack:
        lifted_bed = str(round_dir / f"{st.base}_round{round_num}_detections_original.bed")
        run_lift(round_bed, st.taco_stack, lifted_bed,
                 log_name=f"taco_lift_round{round_num}_detections{tag}")
        st.lifted_beds.append(lifted_bed)
    else:
        st.lifted_beds.append(round_bed)

    orig_elem = _parse_lifted_bed(st.lifted_beds[-1])

    orig_flank = {}
    if st.taco_stack:
        flank_bed = str(round_dir / f"{st.base}_round{round_num}_flank.bed")
        write_flank_bed(extracted_fa, st.current_genome, flank_bed)
        flank_lifted = str(round_dir / f"{st.base}_round{round_num}_flank_original.bed")
        run_lift(flank_bed, st.taco_stack, flank_lifted,
                 log_name=f"taco_lift_round{round_num}_flank{tag}")
        orig_flank = _parse_lifted_bed(flank_lifted)

    return orig_elem, orig_flank


def _reconcile_nodes(st, extracted_fa, orig_elem, orig_flank, round_num,
                     cross_genome, skip_converged, no_node_reuse):
    """Compare this round's detections into the genome's persistent node set.
    """
    st.coords_current = {}

    raw    = _read_fasta_records(extracted_fa)        # orig_name -> (header, seq)
    coords = parse_element_coordinates(extracted_fa)  # orig_name -> (scaf, s, e)
    by_contig: dict = {}
    for orig_name, (header, seq) in raw.items():
        scaffold = coords.get(orig_name, ('', 0, 0))[0]
        if not scaffold:
            continue
        pname   = st.node_name(orig_name, round_num, cross_genome)
        pheader = st.node_header(header, orig_name, pname, cross_genome)
        by_contig.setdefault(scaffold, {})[pname] = (orig_name, pheader, seq)

    cur_order = scaffold_order(st.current_genome)
    if st.scaffold_count is None:
        st.scaffold_count = len(cur_order)
    elif len(cur_order) != st.scaffold_count:
        raise RuntimeError(
            f"{st.base}: scaffold count changed {st.scaffold_count} → "
            f"{len(cur_order)} at round {round_num}; seq_index reuse keys "
            f"would be invalid")

    if skip_converged and st.next_active_file and Path(st.next_active_file).exists():
        with open(st.next_active_file) as af:
            reran_contigs = {l.strip() for l in af if l.strip()}
    else:
        reran_contigs = set(cur_order)

    new_names: set = set()
    if round_num == 1 or not st.taco_stack or no_node_reuse:
        for scaffold in sorted(reran_contigs):
            for old in st.contig_nodes.get(scaffold, set()):
                st.node_seq.pop(old, None)
                st.node_orig_coords.pop(old, None)
            st.contig_nodes[scaffold] = set()
            for pname, (orig_name, pheader, seq) in by_contig.get(scaffold, {}).items():
                st.node_seq[pname] = (pheader, seq)
                st.contig_nodes[scaffold].add(pname)
                st.coords_current[pname] = coords[orig_name]
                oc = orig_elem.get(orig_name)
                if oc is not None:
                    st.node_orig_coords[pname] = oc
                new_names.add(pname)
    else:
        # Incremental reuse path
        n_reused = 0
        for scaffold in sorted(reran_contigs):
            old_by_coord: dict = {}
            for old in sorted(st.contig_nodes.get(scaffold, set())):
                oc = st.node_orig_coords.get(old)
                if oc is not None:
                    old_by_coord.setdefault(oc, []).append(old)
            survivors: set = set()
            for pname, (orig_name, pheader, seq) in by_contig.get(scaffold, {}).items():
                es = orig_elem.get(orig_name)
                fs = orig_flank.get(orig_name)
                affected = False
                if fs is not None:
                    si, fa, fb = fs
                    for ra, rb in st.last_removed_orig.get(si, ()):  # half-open
                        if fa < rb and ra < fb:
                            affected = True
                            break
                match = None
                if not affected and es is not None:
                    cands = old_by_coord.get(es)
                    if cands:
                        for i, cand in enumerate(cands):
                            if st.node_seq.get(cand, (None, None))[1] == seq:
                                match = cand
                                cands.pop(i)
                                break
                if match is not None:
                    st.coords_current[match] = coords[orig_name]
                    survivors.add(match)
                    n_reused += 1
                else:
                    st.node_seq[pname] = (pheader, seq)
                    st.coords_current[pname] = coords[orig_name]
                    if es is not None:
                        st.node_orig_coords[pname] = es
                    new_names.add(pname)
                    survivors.add(pname)
            # drop old nodes neither reused nor re-detected on this contig
            for old in st.contig_nodes.get(scaffold, set()):
                if old not in survivors:
                    st.node_seq.pop(old, None)
                    st.node_orig_coords.pop(old, None)
            st.contig_nodes[scaffold] = survivors
        log.info(f"  [incremental] reused {n_reused} node(s), {len(new_names)} new "
                 f"across {len(reran_contigs)} active contig(s)")

    return new_names

# ---------------------------------------------------------------------------
# One round, all genomes: pooled universe, clustering, classification
# ---------------------------------------------------------------------------

def _build_universe(states, round_dir, round_num, base, new_names):
    """Pool every genome's current nodes into one record set for clustering."""
    from engine.feature_cluster import _load_records

    universe_fa = str(round_dir / f"{base}_round{round_num}_universe.fasta")
    with open(universe_fa, 'w') as out:
        for st in states:
            for name, (header, seq) in st.node_seq.items():
                out.write(header + '\n')
                out.write(seq + '\n')
    records = _load_records(universe_fa)

    Path(universe_fa).unlink(missing_ok=True)
    n_nodes     = sum(len(st.node_seq) for st in states)
    n_collapsed = sum(len(st.collapsed_nodes) for st in states)
    n_stable    = n_nodes - len(new_names) - n_collapsed
    log.info(f"  Universe: {len(records)} nodes "
             f"({len(new_names)} new, {n_collapsed} collapsed-back, "
             f"{n_stable} stable)")
    return records


def _cluster_round(records, states, new_names, reuse_store, frozen_params,
                   fcfg, kmer_lib, wfa_lib, clustering_dir, round_num,
                   no_node_reuse):
    """Cluster the pooled detections, reusing surviving edges.
    """
    from engine.feature_cluster import cluster_records_incremental, _write_gml

    live = set()
    for st in states:
        live.update(st.node_seq)
    reuse_edges = {(a, b): attrs for (a, b), attrs in reuse_store.items()
                   if a in live and b in live}
    first_round = (round_num == 1 or not reuse_edges or no_node_reuse)
    with timed(f"round {round_num}: incremental cluster"):
        consensi_fa, membership_tsv, _all_results, G = cluster_records_incremental(
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
    return consensi_fa, membership_tsv, G, reuse_store


def _classify_round(consensi_fa, membership_tsv, clustering_dir, round_num,
                    args, fcfg):
    """One TEsorter pass over the round's consensi.
    """
    from engine.feature_cluster import write_classified_gml

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
    return classifications, membership, classified_fa


# ---------------------------------------------------------------------------
# One round, one genome: collapse
# ---------------------------------------------------------------------------

def _collapse_genome(st, membership, classifications, round_num, max_rounds,
                     round_dir, keep_genomes, cross_genome):
    """Collapse this genome's LTR-classified detections out of its own sequence.
    """
    added_coords    = []
    collapsed_names: set = set()
    for cons_name, members in membership.items():
        if not _is_ltr_class(classifications.get(cons_name, 'unknown')):
            continue
        for m in members:
            if m in st.coords_current:     # only THIS genome's fresh detections
                added_coords.append(st.coords_current[m])
                collapsed_names.add(m)

    n_events = 0
    st.last_removed_orig = {}
    if round_num < max_rounds and added_coords:
        round_dir.mkdir(parents=True, exist_ok=True)
        tag = st.tool_tag(cross_genome)
        spec_tsv = str(round_dir / f"{st.base}_round{round_num}_collapse.tsv")
        output_prefix = str(round_dir / f"{st.base}_round{round_num}")
        with timed(f"round {round_num}: build collapse spec"):
            n_events, active_scaffolds, n_dropped = _collapse_spec_from_coords(
                added_coords, st.current_genome, spec_tsv,
            )
        if n_dropped:
            log.warn(f"{st.base}: {n_dropped} collapse coordinate(s) named a "
                     f"scaffold absent from {Path(st.current_genome).name} and "
                     f"were dropped — these elements will not be collapsed.")
        st.next_active_file = str(round_dir / f"active_round{round_num + 1}.txt")
        with open(st.next_active_file, 'w') as af:
            for name in sorted(active_scaffolds):
                af.write(name + '\n')

        if n_events > 0:
            # Lift this round's removed intervals to original space for the next round.
            if st.taco_stack:
                removed_lifted = str(
                    round_dir / f"{st.base}_round{round_num}_removed_original.bed")
                run_lift(spec_tsv, st.taco_stack, removed_lifted,
                         log_name=f"taco_lift_round{round_num}_removed{tag}")
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
                    st.last_removed_orig.setdefault(si, []).append((ra, rb))

            with timed(f"round {round_num}: taco collapse"):
                collapsed_fa, taco_file = run_collapse(
                    st.current_genome, spec_tsv, output_prefix,
                    log_name=f"taco_collapse_round{round_num}{tag}",
                )
            st.taco_stack.append(taco_file)
            if PRUNE_GENOMES:
                _prune_genome_snapshot(st.current_genome, keep=keep_genomes)
            st.current_genome = collapsed_fa
            new_size = seq_size(collapsed_fa)
            log.step(f"  Collapsed {n_events:,} interval(s) → genome {new_size:,} bp "
                     f"({100.0 * (1 - new_size / st.original_size):.1f}% removed)")

            # Collapsed detections become frozen nodes: out of the genome,
            # but kept as clustering evidence across future rounds.
            for m in collapsed_names:
                scaf = st.coords_current[m][0]
                st.contig_nodes.get(scaf, set()).discard(m)
                st.collapsed_nodes.add(m)

    if round_num < max_rounds and not added_coords:
        st.converged = True
        tag = f" [{st.base}]" if cross_genome else ''
        log.stage(f"  No new LTR members to collapse — converged.{tag}")

    return n_events, collapsed_names


def _round_stats_row(st, round_num, n_events, n_collapsed, n_ltr_fams, n_non_ltr,
                     t_round, cross_genome):
    genome_size_now = seq_size(st.current_genome)
    row = {
        'round': round_num,
        'detections': st.n_new,
        'families': n_ltr_fams,
        'non_ltr': n_non_ltr,
        'new_ltr_members': n_collapsed,
        'collapse_events': n_events,
        'matched': '-',
        'genome_size_bp': genome_size_now,
        'compression_pct': 100.0 * (1 - genome_size_now / st.original_size),
        'active_pct': st.round_active_pct,
        'skipped_bp': st.round_skipped_bp,
        'time_s': time.perf_counter() - t_round,
    }
    if cross_genome:
        row['genome'] = st.base
    return row


# ---------------------------------------------------------------------------
# Final outputs
# ---------------------------------------------------------------------------

def _finalize(states, output_dir, base, final_classified_src,
              final_membership_src, round_stats, cross_genome):
    """Write all_rounds/ and final/, and print the summary table."""
    all_rounds_dir = Path(output_dir) / 'all_rounds'
    all_rounds_dir.mkdir(parents=True, exist_ok=True)

    for st in states:
        stem = st.base if cross_genome else base
        combined_bed = str(all_rounds_dir / f"{stem}_all_rounds_original.bed")
        with open(combined_bed, 'w') as out:
            out.write("#seq_index\tstart\tend\telement_id\tround\n")
            for rnd, bed in enumerate(st.lifted_beds, start=1):
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
    if not cross_genome:
        log.info(f"  all_rounds/ detections:   "
                 f"{all_rounds_dir / f'{base}_all_rounds_original.bed'}")
    else:
        log.info(f"  all_rounds/ detections:   one BED per genome under "
                 f"{all_rounds_dir}")

    ltr_fa, rm_fa, n_kept, n_dropped = write_final_library(
        final_classified, output_dir, base)
    log.info(f"  {n_dropped} consensi dropped as non-LTR, {n_kept} kept")
    log.stage(f"\n  {ltr_fa} is the final output LTR-RT library "
              f"({n_kept} families)")
    log.stage(f"  {rm_fa}   (RepeatMasker: >NAME#CLASS)") # writing an additional file compatible with RepeatMasker

    _print_incremental_summary(round_stats, cross_genome)
    return final_classified, final_membership


# ---------------------------------------------------------------------------
# Main iterative loop
# ---------------------------------------------------------------------------

def _run_round_loop(states, base, output_dir, args, cross_genome):
    """Iterative collapse over one or more genomes sharing one clustering.
    """
    # libgomp reads OMP_NUM_THREADS exactly once, in its load-time constructor,
    # so this MUST run before anything for parallelism
    os.environ['OMP_NUM_THREADS'] = str(max(1, args.threads))

    from engine.feature_cluster import _ALGO_DIR
    from algorithms.c_libs import KmerLib, WFALib

    kmer_lib = KmerLib(_ALGO_DIR / 'kmer_seed.so')
    wfa_lib  = WFALib(_ALGO_DIR / 'wfa_align.so')
    fcfg     = build_feature_config(args)

    # ── state shared by every genome ────────────────────────────────────────
    reuse_store:   dict = {}   # (name_a, name_b) -> edge attrs (surviving edges)
    frozen_params: dict = {}   # set-dependent seeding params, frozen at round 1
    round_stats:   list = []

    keep_genomes = set()
    for st in states:
        keep_genomes.add(str(st.current_genome))
        keep_genomes.add(str(st.source_path or ''))

    skip_converged = True

    no_node_reuse = bool(os.environ.get('HIREPEAT_NO_NODE_REUSE'))
    if no_node_reuse:
        log.stage("  [HIREPEAT_NO_NODE_REUSE] node + edge reuse DISABLED "
                  "(full universe rebuild + from-scratch re-cluster each round)")

    final_classified_src = None
    final_membership_src = None

    # ── Round 0: tandem compression layers ──────────────────────────────────
    for st in states:
        _round0_tandem(st, st.tandem_layers, round_stats, cross_genome)

    max_rounds = args.rounds

    for round_num in range(1, max_rounds + 1):
        round_dir = output_dir / f"round{round_num}"
        round_dir.mkdir(parents=True, exist_ok=True)

        active = [st for st in states if not st.converged]
        if not active:
            break
        if cross_genome:
            log.stage(f"\n══ Round {round_num}/{max_rounds} ── "
                      f"{len(active)}/{len(states)} genome(s) active")

        # ── per genome: detect, lift, reconcile ─────────────────────────────
        t_round   = {}
        new_names = set()
        for st in active:
            t_round[st.base] = time.perf_counter()
            gdir = st.round_dir(output_dir, round_num, cross_genome)
            extracted_fa = _detect_round(st, round_num, max_rounds, gdir, args,
                                         skip_converged, cross_genome)
            orig_elem, orig_flank = _lift_round_beds(
                st, extracted_fa, round_num, gdir, cross_genome)
            new_names |= _reconcile_nodes(
                st, extracted_fa, orig_elem, orig_flank, round_num,
                cross_genome, skip_converged, no_node_reuse)
        for st in states:
            if st.converged:          # contributes nodes, but nothing collapsible
                st.coords_current = {}
                st.last_removed_orig = {}

        if cross_genome:
            seen: dict = {}
            for st in states:
                stray = set(st.coords_current) - set(st.node_seq)
                if stray:
                    raise RuntimeError(
                        f"{st.base}: {len(stray)} collapse-eligible name(s) are "
                        f"not nodes of this genome, e.g. {sorted(stray)[:3]}")
                for name in st.coords_current:
                    if name in seen:
                        raise RuntimeError(
                            f"element {name} is collapse-eligible in both "
                            f"{seen[name]} and {st.base}")
                    seen[name] = st.base

        # ── pooled universe → one clustering → one classification ───────────
        records = _build_universe(states, round_dir, round_num, base, new_names)
        clustering_dir = round_dir / 'clustering'
        consensi_fa, membership_tsv, _G, reuse_store = _cluster_round(
            records, states, new_names, reuse_store, frozen_params,
            fcfg, kmer_lib, wfa_lib, clustering_dir, round_num, no_node_reuse,
        )
        classifications, membership, classified_fa = _classify_round(
            consensi_fa, membership_tsv, clustering_dir, round_num, args, fcfg)

        final_classified_src = classified_fa if Path(classified_fa).exists() else consensi_fa
        final_membership_src = membership_tsv

        n_ltr_fams = sum(1 for cls in classifications.values() if _is_ltr_class(cls))
        n_non_ltr  = len(classifications) - n_ltr_fams
        log.step(f"  Round {round_num}: {len(membership)} families "
                 f"({n_ltr_fams} LTR, {n_non_ltr} non-LTR)")

        # ── per genome: collapse its own members of the LTR families ────────
        for st in active:
            gdir = st.round_dir(output_dir, round_num, cross_genome)
            n_events, collapsed_names = _collapse_genome(
                st, membership, classifications, round_num, max_rounds,
                gdir, keep_genomes, cross_genome)
            round_stats.append(_round_stats_row(
                st, round_num, n_events, len(collapsed_names),
                n_ltr_fams, n_non_ltr, t_round[st.base], cross_genome))

        if all(st.converged for st in states):
            break

    final_classified, final_membership = _finalize(
        states, output_dir, base, final_classified_src, final_membership_src,
        round_stats, cross_genome)

    taco_stacks = [st.taco_stack for st in states]
    return ([final_classified], [final_membership],
            taco_stacks[0] if len(taco_stacks) == 1 else taco_stacks,
            round_stats)


def run_rounds(genome_for_ltr, base, output_dir, round1_extracted,
               original_genome_size, args, tandem_layers=None):
    """Exact incremental TE collapse via a persistent element similarity graph.

    Single-genome: one GenomeState through the shared loop.
    """
    st = GenomeState(
        source_path=str(getattr(args, 'genome', '') or ''),
        base=base,
        original_size=original_genome_size,
        current_genome=str(genome_for_ltr),
        round1_extracted=str(round1_extracted),
        tandem_layers=list(tandem_layers or []),
    )
    return _run_round_loop([st], base, output_dir, args, cross_genome=False)


def make_genome_state(genome_path, base, original_size, genome_for_ltr,
                      round1_extracted, tandem_layers):
    """One genome's state for the cross-genome collapse loop."""
    return GenomeState(
        source_path=str(genome_path),
        base=base,
        original_size=original_size,
        current_genome=str(genome_for_ltr),
        round1_extracted=str(round1_extracted),
        tandem_layers=list(tandem_layers or []),
    )


def run_rounds_crossgenome(states, base, output_dir, args):
    """Iterative collapse across several genomes sharing one clustering.
    """
    log.stage(f"\nIterative collapse over {len(states)} genome(s), "
              f"up to {args.rounds} round(s), clustered together")
    return _run_round_loop(states, base, output_dir, args, cross_genome=True)


def _print_incremental_summary(round_stats, cross_genome=False):
    """Per-round summary table for the incremental loop."""
    gcol = (f" {'Genome':<16}" if cross_genome else '')
    gsep = (f" {'-'*16}" if cross_genome else '')
    log.stage(f"\n{'=' * 80}")
    log.stage("  Incremental collapse summary")
    log.stage(f"{'=' * 80}")
    log.stage(f"  {'Round':<9}{gcol} {'Detect':>7} {'Matched':>8} {'LTRfam':>7} "
              f"{'non-LTR':>8} {'NewMem':>7} {'Collapsed':>10} {'Genome_Mbp':>11} "
              f"{'Active%':>8} {'Time(s)':>8}")
    log.stage(f"  {'-'*9}{gsep} {'-'*7} {'-'*8} {'-'*7} {'-'*8} {'-'*7} {'-'*10} "
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
        gval   = (f" {str(s.get('genome', '-')):<16}" if cross_genome else '')
        log.stage(f"  {str(s['round']):<9}{gval} {str(det):>7} {str(mat):>8} {str(ltr_):>7} "
                  f"{str(nlt):>8} {str(newm):>7} {str(coll):>10} {gbp_s:>11} "
                  f"{apct_s:>8} {t_s:>8.1f}")
    log.stage('')
