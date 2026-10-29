"""Genome-path naming helpers and output-directory setup."""

import re
from pathlib import Path

from utilities import log


def genome_label(path):
    """Sample ID with all extensions stripped, e.g. 'idAnoFuneDA-414_04'."""
    name = Path(path).name
    while '.' in name:
        name = name.rsplit('.', 1)[0]
    return name


def strip_genome_name(path):
    """Sample ID plus its .pri/.alt haplotype tag, if the path carries one."""
    name = Path(path).name
    hap = next((f".{tag}" for tag in ("pri", "alt") if f".{tag}." in name), "")
    return f"{genome_label(path)}{hap}"


def gml_safe_label(label):
    return re.sub(r'[-_.]', '', label)


def set_folder_structure(path, suffix="_ltr_output", parent=None):
    """Create <stem><suffix>/ in *parent* (default CWD).  Returns (stem, Path)."""
    name = strip_genome_name(path)
    output_dir = Path(parent or '.') / f"{name}{suffix}"
    output_dir.mkdir(parents=True, exist_ok=True)
    log.stage(f"Output directory: {output_dir}")
    return name, output_dir
