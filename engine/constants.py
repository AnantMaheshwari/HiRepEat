"""Tuning constants for the HiRepEat pipeline.
"""

# ---------------------------------------------------------------------------
# Clustering: the Wicker 80/80/80 rule
# ---------------------------------------------------------------------------

MIN_IDENTITY = 0.80
MIN_COVERAGE = 0.80

# Smallest connected component kept as a family.
MIN_CLUSTER_SIZE = 3

# discard k-mers occurring in more than this many distinct elements (LTRs)
MAX_BUCKET = 50

# Stricter for internal regions
MAX_BUCKET_INTERNAL = 25

# Flanking-sequence ortholog/tandem filter.
FLANK_FILTER = True

# Shortest flank worth comparing.
MIN_FLANK_LEN = 50


# ---------------------------------------------------------------------------
# k-mer seeding and WFA
# ---------------------------------------------------------------------------

KMER_SIZE = 12

# Minimum shared k-mers before a pair becomes a WFA candidate.
MIN_SEED_SHARED = 3

# WFA edit budget as a fraction of max(la, lb).
MAX_EDIT_FRAC = 0.25

# Shortest feature regions worth including in a cmparison.
MIN_LTR_LEN      = 80
MIN_INTERNAL_LEN = 300


# ---------------------------------------------------------------------------
# Consensus building
# ---------------------------------------------------------------------------

# Per-column base frequency needed to call a consensus base; below it the column
# becomes N, or an IUPAC code when AMBIGUITY_THRESHOLD is enabled.
CONS_THRESHOLD = 0.4

# Columns occupied in fewer than this fraction of members are dropped.
SATURATION_THRESHOLD = 0.6

# Lower bar for emitting an IUPAC ambiguity code when no single base reaches
# CONS_THRESHOLD.  0 disables it: the column just becomes N.
AMBIGUITY_THRESHOLD = 0.0

# WFA budget for one progressive profile merge, as a fraction of the longer
# child consensus.
MERGE_MAX_EDIT_FRAC = 0.20

# Identity below which a progressive merge is refused and the subtree splits.
MIN_MERGE_IDENTITY = 0.80

# Split a family into subfamilies by cutting the UPGMA guide tree: merge steps
# at or below CUT_THRESHOLD stay in one subfamily.  0.20 ~ 80% identity.
SPLIT_SUBFAMILIES = True
CUT_THRESHOLD     = 0.20

# Outgroup tie-breaking: at each internal node, use the sibling subtree's consensus to
# break near-ties in ambiguous columns
USE_OUTGROUP = True

# Same as max_bucket but for consensus building.
CONSENSUS_ANCHOR_MAX_BUCKET = 50


# ---------------------------------------------------------------------------
# Tandem compression and iterative collapse
# ---------------------------------------------------------------------------

# FasTAN iteration's -r restriction, must be >= FasTAN's DIAG_MAX (30000)
# Allows iteratve tandem search to restrict to just the parts of contig that change
FASTAN_SCAN_FLANK = 32768

# Delete intermediate .1seq files, everything is still reconstructible but for most 
# cases the intermediate .1seqs are not useful, only the tandem-deleted and final .1seq
# need to be stored alongside the LTR-RT library (space optimization for running at a larger scale)
PRUNE_GENOMES = True


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

# TEsorter HMM database.  'rexdb' is REXdb viridiplantae v4.0 + metazoa v3.1.
TESORTER_DB = 'rexdb'


def tesorter_gydb() -> bool:
    """Whether to run a second TEsorter pass under gydb.
    """
    return TESORTER_DB != 'gydb'
