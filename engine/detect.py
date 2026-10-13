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
    genome_label, set_folder_structure, strip_genome_name,
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

    return args


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
        check_resume, cluster_and_classify, combine_genome_fastas,
        detect_and_extract, guard_tandem_only_outputs, run_rounds,
        run_tandem_only,
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
        output_dir = list_path.parent / f"{base}_ltr_output"
        output_dir.mkdir(exist_ok=True)
        setup_logging(output_dir / f"{base}_pipeline.log")

        if args.rounds == 0:
            stems = [Path(g).stem for g in genome_paths]
            dupes = sorted({s for s in stems if stems.count(s) > 1})
            if dupes:
                raise SystemExit(
                    f"--rounds 0: genomes sharing a tandem/ prefix would "
                    f"overwrite each other: {', '.join(dupes)}")
            for gpath in genome_paths:
                guard_tandem_only_outputs(
                    output_dir, gpath,
                    output_dir / f"{strip_genome_name(gpath)}_ltr_sequences.fasta")

        write_run_params(output_dir / 'run_params.yaml', args, output_dir, base,
                         tool_versions=tool_versions)

        t_pipeline = time.perf_counter()
        log.stage(f"Genome list: {args.genome_list}  ({len(genome_paths)} genomes)")
        log.stage(f"Output directory: {output_dir}")

        if args.rounds == 0:
            for gpath in genome_paths:
                log.stage(f"\n── Genome: {gpath}  (label: {genome_label(gpath)})")
                run_tandem_only(gpath, output_dir, args.tandem_rounds,
                                prune_genomes=constants.PRUNE_GENOMES)
        else:
            extracted_fastas = []   # (Path, label)
            for gpath in genome_paths:
                label = genome_label(gpath)
                log.stage(f"\n── Genome: {gpath}  (label: {label})")
                efasta, _, _ = detect_and_extract(
                    gpath, strip_genome_name(gpath), output_dir, args.tandem_rounds,
                    prune_genomes=constants.PRUNE_GENOMES,
                )
                extracted_fastas.append((efasta, label))

            combined_fasta = output_dir / f"{base}_combined_ltr_sequences.fasta"
            log.step(f"Combining {len(extracted_fastas)} genome FASTAs → {combined_fasta}")
            combine_genome_fastas(extracted_fastas, combined_fasta)

            consensus_fasta, classified_fasta = cluster_and_classify(
                str(combined_fasta), output_dir / 'clustering', args,
            )
            log.stage(f"Consensus library saved to {consensus_fasta}")
            log.stage(f"Classified consensus library saved to {classified_fasta}")

            # final/ here too
            from utilities.final_library import write_final_library
            ltr_fa, rm_fa, n_kept, n_dropped = write_final_library(
                classified_fasta, output_dir, base)
            log.info(f"  {n_dropped} consensi dropped as non-LTR, {n_kept} kept")
            log.stage(f"\n  {ltr_fa} is the final output LTR-RT library "
                      f"({n_kept} families)")
            log.stage(f"  {rm_fa}   (RepeatMasker: >NAME#CLASS)")

    else:
        # ── Single-genome path ────────────────────────────────────────────────
        base, output_dir = set_folder_structure(args.genome, "_ltr_output")
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
