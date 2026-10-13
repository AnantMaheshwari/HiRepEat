"""
Classify LTR consensus sequences using TEsorter (rexdb HMM database).

Runs TEsorter — from PATH, or out of the HiRepEat container image — and parses
the .cls.tsv output to produce a classified FASTA with class/domain tags
appended to each header.
"""

import time
from collections import Counter
from pathlib import Path
from typing import List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter

from utilities import log
from utilities.external import tesorter_cmd, tool_env
from utilities.final_library import is_ltr
from utilities.timing import timed


_UNCLASSIFIED = ('Unknown', 'none', '?', 'no', 'unknown')

# GyDB2 and REXdb different terminology for same superfamily
_SUPERFAMILY_ALIASES = {
    'Pao': 'Bel-Pao',
    'Retroviridae': 'Retrovirus',
}


def _as_db_list(db) -> List[str]:
    """Accept either a single database name or a sequence of them."""
    if isinstance(db, str):
        return [db]
    return list(db)


def run_tesorter(
    input_fasta: str,
    output_dir: str,
    db: str = "rexdb",
    threads: int = 4,
    log_name: str = "tesorter",
    prefix: str = None,
    tmp_dir: str = None,
) -> str:
    """Run TEsorter on input FASTA.  Returns the path to the .cls.tsv output."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if prefix is None:
        prefix = str(out_dir / Path(input_fasta).name)
    if tmp_dir is None:
        tmp_dir = str(out_dir / 'tesorter_tmp')

    cmd = tesorter_cmd(
        [input_fasta,
         '-db', db,
         '-p', threads,
         '-pre', prefix,
         '-tmp', tmp_dir,
         '-fw'],
        bind_paths=(input_fasta, out_dir),
    )
    log.info(f"  Running: TEsorter -db {db} -p {threads}")
    log.run_tool(cmd, log_name, env=tool_env())

    cls_tsv = f"{prefix}.cls.tsv"
    if not Path(cls_tsv).exists():
        raise FileNotFoundError(
            f"TEsorter did not produce expected output: {cls_tsv}")
    return cls_tsv


def parse_tesorter_tsv(cls_tsv: str) -> dict:
    """Parse TEsorter .cls.tsv into {seq_id: (class_label, domains_str, strand, complete)}.

    Maps TEsorter Order/Superfamily to class labels:
      LTR/Gypsy, LTR/Copia, LTR/Bel-Pao, LTR/ERV, LTR/unknown, Unknown

    """
    results = {}
    with open(cls_tsv) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            te_id = parts[0]
            order = parts[1]          # e.g. "LTR"
            superfamily = parts[2]    # e.g. "Gypsy", "Copia", "Bel-Pao"
            superfamily = _SUPERFAMILY_ALIASES.get(superfamily, superfamily)
            clade = parts[3]          # e.g. "unknown", "Athila", etc.
            complete = parts[4]       # "yes" or "no"
            strand = parts[5]         # "+" or "-"
            domains = parts[6]        # e.g. "GAG|Ty1_copia RT|Ty1_copia ..."

            # Build class label
            if order == "LTR" and superfamily != "unknown":
                class_label = f"LTR/{superfamily}"
            elif order != "Unknown":
                class_label = f"{order}/unknown"
            else:
                class_label = "Unknown"

            results[te_id] = (class_label, domains, strand, complete, clade)

    return results


def _merge_calls(seq_id: str, per_db: List[Tuple[str, dict]]):
    """Pick one call for `seq_id` across the databases, in the order given.
    """
    primary_db, primary_calls = per_db[0]
    primary = primary_calls.get(seq_id)
    if primary is not None and is_ltr(primary[0]):
        return primary, primary_db

    for db, calls in per_db[1:]:
        call = calls.get(seq_id)
        if call is not None and is_ltr(call[0]):
            return call, db

    if primary is not None:
        return primary, primary_db
    for db, calls in per_db[1:]:
        call = calls.get(seq_id)
        if call is not None:
            return call, db
    return _UNCLASSIFIED, None


def classify_consensi_tesorter(
    input_fasta: str,
    output_fasta: str,
    db="rexdb",
    threads: int = 4,
    log_name: str = "tesorter",
) -> Tuple[List[str], List[str]]:
    """Classify consensus sequences using TEsorter.

    Parameters
    ----------
    input_fasta : str
        Input FASTA with consensus sequences.
    output_fasta : str
        Output FASTA with class=, domains=, strand=, clade= tags appended.
    db : str or sequence of str
        TEsorter HMM database name(s)
    threads : int
        Number of processors for TEsorter.

    Returns
    -------
    predicted_names : list of str
    domain_strings : list of str
    """
    t_total = time.perf_counter()
    dbs = _as_db_list(db)

    # Load sequences
    log.info(f"  Loading sequences from {input_fasta} ...")
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        log.info("  No sequences to classify.")
        with open(output_fasta, "w"):
            pass
        return [], []

    log.step("  TEsorter classification ...")
    log.info(f"  {len(records)} consensi, -db {' + '.join(dbs)}")

    out_dir = Path(output_fasta).parent


    per_db: List[Tuple[str, dict]] = []
    for i, db_name in enumerate(dbs):
        suffix = '' if i == 0 else f'.{db_name}'
        with timed(f"TEsorter ({db_name})"):
            cls_tsv = run_tesorter(
                input_fasta, str(out_dir), db=db_name, threads=threads,
                log_name=log_name if i == 0 else f"{log_name}_{db_name}",
                prefix=str(out_dir / (Path(input_fasta).name + suffix)),
                tmp_dir=str(out_dir / f"tesorter_tmp{suffix}"),
            )
        calls = parse_tesorter_tsv(cls_tsv)
        per_db.append((db_name, calls))
        n_ltr = sum(1 for c in calls.values() if is_ltr(c[0]))
        log.info(f"  {len(calls)}/{len(records)} consensi classified by "
                 f"TEsorter -db {db_name} ({n_ltr} as LTR)")

    # Build classified records
    predicted_names = []
    domain_strings = []
    classified_records = []
    won_by = Counter()

    for rec in records:
        (class_label, domains, strand, complete, clade), src_db = _merge_calls(
            rec.id, per_db)
        if src_db is not None:
            won_by[src_db] += 1

        predicted_names.append(class_label)
        domain_strings.append(domains)

        new_desc = (
            f"{rec.description} class={class_label} "
            f"domains={domains} strand={strand} "
            f"clade={clade} complete={complete}"
        )
        if len(dbs) > 1:
            new_desc += f" db={src_db or 'none'}"
        classified_records.append(
            SeqRecord(rec.seq, id=rec.id, description=new_desc))

    n_classified = sum(won_by.values())
    log.info(f"  {n_classified}/{len(records)} consensi classified by TEsorter")

    counts = Counter(predicted_names)
    log.info("  Classification results: "
             + ', '.join(f"{cls} {n}" for cls, n in counts.most_common()))

    if len(dbs) > 1:
        primary = per_db[0][1]
        n_primary_ltr = sum(1 for c in primary.values() if is_ltr(c[0]))
        n_union_ltr = sum(1 for c in predicted_names if is_ltr(c))
        log.info(f"  Union of {len(dbs)} databases: {n_union_ltr} LTR calls "
                 f"({n_union_ltr - n_primary_ltr:+d} vs -db {dbs[0]} alone)")
        log.info("  Calls kept per database: "
                 + ', '.join(f"{d} {won_by.get(d, 0)}" for d in dbs))

    # Write output FASTA
    with open(output_fasta, 'w') as _h:
        FastaWriter(_h, wrap=None).write_file(classified_records)
    log.info(f"  Wrote {len(classified_records)} classified sequences to "
             f"{output_fasta}")

    log.info(f"  [timing] classify_consensi_tesorter total: "
             f"{time.perf_counter() - t_total:.2f}s")

    return predicted_names, domain_strings
