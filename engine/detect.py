"""HiRepEat pipeline wrapper + argument parsing, calls into pipeline.py
"""

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import argparse
import time
from pathlib import Path

from engine import constants
from seq.info import seq_size
from utilities import log
from utilities.external import TESORTER_DBS, preflight, set_image
from utilities.log import setup_logging
from utilities.path_operations import (
    gml_safe_label, set_folder_structure, strip_genome_name,
)
from utilities.provenance import read_run_params, write_run_params


def parse_args():
    parser = argparse.ArgumentParser(
        description='Hierarchical Repeat Eater: LTR-RT detection, clustering, classification, and iterative collapse.'
    )
    parser.add_argument('--genome',
                        help='Path to primary genome FASTA (can be gzipped). '
                             'Required unless --genome-list is given.')
    parser.add_argument('--genome-list',
                        dest='genome_list',
                        help='Path to a text file listing genome FASTA paths (one per line). '
                             'Each genome is processed independently through detection and '
                             'extraction, then all extracted sequences are clustered together '
                             'into one shared library (cross-genome mode).')
    parser.add_argument('--output', default=None,
                        help='Directory to write the <genome>_ltr_output/ folder into '
                             '(created if missing). Default: the current directory, '
                             "or the --genome-list file's directory in cross-genome mode.")
    parser.add_argument('--tandem-rounds', type=int, default=10,
                        dest='tandem_rounds',
                        help='Max number of iterative FasTAN tandem compression rounds. '
                             'Each round runs FasTAN then taco compress on the result. '
                             'Stops early when a round finds no new tandems. '
                             '0 disables tandem compression entirely (default: 10)')
    parser.add_argument('--rounds', type=int, default=10,
                        help='Number of iterative collapse rounds. Round 1 is the initial '
                             'FastLTR detection. Rounds 2+ collapse detected LTR elements '
                             'from the genome and re-run FastLTR to find elements that were '
                             'previously nested. '
                             '0 runs tandem compression only, then stops. '
                             '(default: 10)')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads for WFA and consensus parallelism (default: 8)')
    parser.add_argument('--debug', action='store_true',
                        help='Write clustering diagnostics: the similarity graphs '
                             '(graph visualization) and the '
                             'per-family WFA-based trees.')
    parser.add_argument('--keep-alignments', action='store_true',
                        dest='keep_alignments',
                        help='Also write the alignment of the cluster members '
                             'behind each consensus, as aligned FASTA, into '
                             '<clustering>/alignments/')
    parser.add_argument('--image', dest='image', default=None,
                        help='Path to the HiRepEat container image (.sif) to '
                             'run TEsorter from. Only needed when TEsorter is '
                             'not on PATH; ignored when already running inside '
                             'the image. Overrides $HIREPEAT_IMAGE.')
    args = parser.parse_args()

    if not args.genome and not args.genome_list:
        parser.error('one of --genome or --genome-list is required')
    if args.genome_list and args.genome:
        parser.error('--genome-list is mutually exclusive with --genome')
    if args.tandem_rounds < 0:
        parser.error('--tandem-rounds must be >= 0')
    if args.rounds < 0:
        parser.error('--rounds must be >= 0')
    if args.rounds == 0 and args.tandem_rounds == 0:
        parser.error('--rounds 0 with --tandem-rounds 0 would do nothing')
    if constants.TESORTER_DB not in TESORTER_DBS:
        parser.error(
            f"engine/constants.py: TESORTER_DB = {constants.TESORTER_DB!r} is "
            f"not a TEsorter database name.\nValid names: "
            + ', '.join(TESORTER_DBS))

    args.flank_filter = constants.FLANK_FILTER and not args.genome_list

    args.genome_paths = []
    if args.genome_list:
        args.genome_paths = [
            l.strip() for l in open(args.genome_list)
            if l.strip() and not l.strip().startswith('#')
        ]
        if not args.genome_paths:
            parser.error(f'--genome-list file is empty: {args.genome_list}')
        _check_unique_genome_names(args.genome_paths, parser)

    return args

# need a few of extra checks for iterative collapse (TODO: Anant refactor all of these into asserts)
def _check_unique_genome_names(paths, parser):
    """Refuse a genome list whose members would share a namespace (iterative multi-genome collapse checks).
    """
    checks = (
        (strip_genome_name, 'sample name'),
        (lambda p: gml_safe_label(strip_genome_name(p)), "GML-safe 'genome=' tag"),
        (lambda p: Path(p).stem, 'tandem-layer file prefix'),
    )
    for key, what in checks:
        seen = {}
        for p in paths:
            seen.setdefault(key(p), []).append(p)
        clashes = {k: v for k, v in seen.items() if len(v) > 1}
        if clashes:
            detail = '\n'.join(f"    {k}: " + ', '.join(v)
                               for k, v in sorted(clashes.items()))
            parser.error(
                f"--genome-list: these genomes share a {what}, so their "
                f"elements, outputs or tandem layers would collide:\n{detail}")


def main():
    from engine.pipeline import ResumeError
    from utilities.external import ExternalToolError
    try:
        _run()
    except (ExternalToolError, ResumeError) as exc:
        log.error(str(exc))
        raise SystemExit(1)


def _run():
    args = parse_args()

    set_image(args.image)
    tool_versions = None if args.rounds == 0 else preflight()

    from engine.pipeline import (
        check_genome_list_resume, check_resume, detect_and_extract,
        guard_tandem_only_outputs, make_genome_state, run_rounds,
        run_rounds_crossgenome, run_tandem_only,
    )

    if args.rounds == 0:
        log.info("--rounds 0: --threads, --debug and --keep-alignments have "
                 "no effect in tandem-only mode.")
    if args.genome_list:
        # ── Cross-genome path ─────────────────────────────────────────────────
        # Detection runs per genome; the extracted sequences are then clustered
        # together as one set. 
        genome_paths = args.genome_paths
        list_path  = Path(args.genome_list)
        base       = list_path.stem
        parent     = Path(args.output) if args.output else list_path.parent
        output_dir = parent / f"{base}_ltr_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        setup_logging(output_dir / f"{base}_pipeline.log")

        check_genome_list_resume(read_run_params(output_dir / 'run_params.yaml'),
                                 genome_paths, output_dir)

        if args.rounds == 0:
            for gpath in genome_paths:
                gbase = strip_genome_name(gpath)
                guard_tandem_only_outputs(
                    output_dir, gpath,
                    [output_dir / f"{gbase}_ltr_sequences.fasta",
                     output_dir / 'round1' / gbase / f"{gbase}_ltr_sequences.fasta"])

        write_run_params(output_dir / 'run_params.yaml', args, output_dir, base,
                         tool_versions=tool_versions)

        t_pipeline = time.perf_counter()
        log.stage(f"Genome list: {args.genome_list}  ({len(genome_paths)} genomes)")
        log.stage(f"Output directory: {output_dir}")

        if args.rounds == 0:
            for gpath in genome_paths:
                log.stage(f"\n── Genome: {gpath}  (label: {strip_genome_name(gpath)})")
                run_tandem_only(gpath, output_dir, args.tandem_rounds,
                                prune_genomes=constants.PRUNE_GENOMES)
        else:
            if (output_dir / 'clustering').is_dir():
                log.warn(f"{output_dir / 'clustering'} is from a previous "
                         f"single-pass run, will be over-written")
            states = []
            for gpath in genome_paths:
                gbase = strip_genome_name(gpath)
                log.stage(f"\n── Genome: {gpath}  (label: {gbase})")
                legacy = output_dir / f"{gbase}_ltr_sequences.fasta"
                detect_dir = output_dir if legacy.exists() else output_dir / 'round1' / gbase
                if legacy.exists():
                    log.step(f"         round-1 detections at {legacy}")
                genome_bp = seq_size(gpath)
                efasta, genome_for_ltr, tandem_layers = detect_and_extract(
                    gpath, gbase, output_dir, args.tandem_rounds,
                    prune_genomes=constants.PRUNE_GENOMES,
                    detect_dir=detect_dir,
                )
                states.append(make_genome_state(
                    gpath, gbase, genome_bp, genome_for_ltr, efasta, tandem_layers))

            run_rounds_crossgenome(states, base, output_dir, args)

    else:
        # ── Single-genome path ────────────────────────────────────────────────
        base, output_dir = set_folder_structure(args.genome, "_ltr_output",
                                                parent=args.output)
        setup_logging(output_dir / f"{base}_pipeline.log")
        check_resume(read_run_params(output_dir / 'run_params.yaml'),
                     args.genome, base, output_dir)

        if args.rounds == 0:
            guard_tandem_only_outputs(
                output_dir, args.genome,
                output_dir / 'round1' / f"{base}_ltr_sequences.fasta")

        original_genome_size = seq_size(args.genome)
        write_run_params(output_dir / 'run_params.yaml', args, output_dir, base,
                         genome_bp=original_genome_size,
                         tool_versions=tool_versions)

        t_pipeline = time.perf_counter()

        log.stage(f"\n── Genome: {args.genome} ({original_genome_size:,} bp)")

        if args.rounds == 0:
            run_tandem_only(args.genome, output_dir, args.tandem_rounds,
                            prune_genomes=constants.PRUNE_GENOMES)
        else:
            extracted_fasta, genome_for_ltr, tandem_layers = detect_and_extract(
                args.genome, base, output_dir, args.tandem_rounds,
                prune_genomes=constants.PRUNE_GENOMES,
                detect_dir=output_dir / 'round1',
            )

            run_rounds(
                genome_for_ltr=genome_for_ltr,
                base=base,
                output_dir=output_dir,
                round1_extracted=str(extracted_fasta),
                original_genome_size=original_genome_size,
                args=args,
                tandem_layers=tandem_layers,
            )

    log.stage(f"Done — total pipeline: {time.perf_counter() - t_pipeline:.2f}s")


if __name__ == '__main__':
    main()
