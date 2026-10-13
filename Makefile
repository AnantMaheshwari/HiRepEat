.PHONY: all external detectors alntools onecode libs seqtools extract_ltrs aln_info sif clean

all: external libs seqtools

# The image is built by CI (.github/workflows/container.yml) and published to
# GHCR; this only pulls it down as a SIF.  Pin by digest for runs you intend to
# keep:  make sif IMAGE=ghcr.io/anantmaheshwari/hirepeat@sha256:<digest>
IMAGE ?= ghcr.io/anantmaheshwari/hirepeat:latest
SIF   ?= HiRepEat.sif

sif:
	singularity pull $(SIF) docker://$(IMAGE)

external: detectors alntools onecode

# Single unified source tree builds FasTAN, FastLTR and FastTIR.
detectors:
	$(MAKE) -C external/FASTAN

alntools:
	$(MAKE) -C external/alntools

onecode:
	$(MAKE) -C external/ONEcode
	$(MAKE) -C external/ONEcode/SEQUENCE_UTILITIES

libs:
	bash algorithms/build_feature_align.sh

# seq/ helper tools — both link against the SEQUENCE_UTILITIES .o files, so they
# depend on `onecode` having been built first.
seqtools: extract_ltrs aln_info

extract_ltrs: onecode
	bash seq/build_extract_ltrs.sh

aln_info: onecode
	bash seq/build_aln_info.sh

clean:
	rm -f algorithms/kmer_seed.so algorithms/wfa_align.so
	rm -f seq/extract_ltrs seq/aln_info
	-$(MAKE) -C external/FASTAN clean
	-$(MAKE) -C external/alntools clean
	-$(MAKE) -C external/ONEcode clean
	-$(MAKE) -C external/ONEcode/SEQUENCE_UTILITIES clean
