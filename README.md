<table>
<tr>
<td width="120">
<img src="images/mosquito_eating_repeats.png" width="120" alt="HiRepEat logo">
</td>
<td>

# <strong>Hi</strong>erarchical<strong>Rep</strong>eat<strong>Eat</strong>er

</td>
</tr>
</table>

<table>
<tr>
<td>

A framework for **de novo discovery, clustering, and classification of LTR retrotransposons**, building on Myers' FasTAN and FastLTR (https://github.com/thegenemyers/FASTAN)

Given a genome (FASTA/.1seq), detects repeat elements, cluster, builds a consensus sequence for each, and classifies them, producing
a TE consensus library suitable for annotation, alongside intermediate outputs including clusters and family phylogenies.

</td>
<td width="220">
<img src="images/cross_genome_high_res.png" width="220" alt="HiRepEat framework overview">
</td>
</tr>
</table>

---

<p align="center">
<img src="images/methods.png" alt="HiRepEat method overview">
</p>

---

```
   GENOME (.fasta/.1seq)
      │
  ┌───▼───────────────────────────────────────────────┐
  │ 1. TANDEM COMPRESSION      FasTAN + taco          │  
  │    compress tandem arrays to one copy each        │  --tandem-rounds N
  └───┬───────────────────────────────────────────────┘
      │  TACO tandem-deleted genome (+ SALSA coordinate maps)
  ┌───▼───────────────────────────────────────────────┐
  │ 2. LTR DETECTION           FastLTR                │
  │    find LTR-retrotransposon copies.               │  (--rounds N)
  └───┬───────────────────────────────────────────────┘
      │  extracted element sequences (5'LTR / internal / 3'LTR / flanks)
  ┌───▼───────────────────────────────────────────────┐
  │ 3. CLUSTERING      k-mer seed + WFA               │
  │ families → subfamily split → consensus sequences  │
  └───┬───────────────────────────────────────────────┘
      │  consensi.fasta
  ┌───▼───────────────────────────────────────────────┐
  │ 4. CLASSIFICATION          TEsorter               │
  └───────────────────────────────────────────────────┘
             │
        classified TE consensus library (.fasta)
```


### 1. Tandem compression

Iterates:
- **FasTAN** tandem repeat detector
- **taco** keeps just **one copy** of each tandem array, and **salsa**
  writes a reversible coordinate map (`.1taco`) so any coordinate can be lifted back
  to the original genome.

`--tandem-rounds N` repeats (FasTAN → taco/salsa) iteratively, in order to delete nested tandem structures.

`--rounds 0` tandem compression only, no LTR detection, clustering, classification, collapse etc.

### 2. LTR detection

```
   ...===LTR===[ internal body ]===LTR===...
```

With `--rounds N`, confidently classified LTR families are collapsed out of the
genome (via **taco/salsa**) and detection is re-run, exposing elements that were previously nested, allowing for their clustering.

### 3. Clustering & consensus building

Similarity graph is built with LTR-RT instances as nodes, edges constructed if LTRs wavefront align to 80% identity, and internal sequence also passes 80/80/80 rule.
**Families** are the connected components of that graph (with >=3 elements).

**Consensus** is built per family by **progressive, UPGMA-style profile merging with outgroup correction** to handle transitivity issues due to the 80/80/80 rule.

**Classification** is done using TEsorter to confirm presence of LTR-RT protein domains (https://github.com/zhangrengang/TEsorter)

## Install / build

### Singularity

**Installation:**

```bash
singularity pull HiRepEat.sif docker://ghcr.io/anantmaheshwari/hirepeat:latest
```

**Execution:**
`singularity exec HiRepEat.sif HiRepEat --genome genome.fa --rounds 10`

The image carries TEsorter, hmmer and BLAST.

### From source

```bash
# 1. Python >= 3.11
python3.11 -m venv env
source env/bin/activate
pip install -r requirements.txt
make            
```

`make` builds everything except TEsorter. Supply that with either the image
above (`export HIREPEAT_IMAGE=/path/to/HiRepEat.sif`), a TEsorter on `PATH`, or
`export HIREPEAT_TOOL_PATH=<tesorter bin>:<blast bin>:<hmmer bin>` for an
existing install.

---

## Usage

### From a container image

```bash
singularity exec HiRepEat.sif HiRepEat --genome genome.fa --threads 8 --rounds 10
```

### From source 

```bash
source env/bin/activate
export PYTHONPATH=$(pwd)
export HIREPEAT_IMAGE=/path/to/HiRepEat.sif    # only if TEsorter is not on PATH

# Simplest run: default rounds: 10
python engine/detect.py --genome genome.fa --threads 8 --rounds 10

# Many genomes clustered together into one shared library (cross-genome mode): default rounds: 1
python engine/detect.py --genome-list genomes.txt --threads 8

# Write results to runs/genome_ltr_output/ instead of the current directory
python engine/detect.py --genome genome.fa --output runs
```

### Outputs
Results go in `<genome>_ltr_output/`, created in the current directory (or next to the list file for `--genome-list`). Use `--output DIR` to specify create desired output location as `DIR` instead.

```final/*_all_rounds_classified.fasta``` contains LTR-RT consensus sequences found, and ```final/*_all_rounds_classified_rm.fasta``` is the same file pre-formatted to give as input to RepeatMasker (https://github.com/Dfam-consortium/RepeatMasker)

## License
MIT, see [LICENSE](LICENSE).

## Citation

If you use HiRepEat in your research, please cite:

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.10.14.759201-b31b1b.svg)](https://doi.org/10.64898/2026.10.14.759201)
> Maheshwari A, Myers EW, Sierra P, Vertebrate
Genomes Project Consortium Phase I, Lawniczak MKN, Durbin R. Efficient and precise discovery of LTR-retrotransposon families using HiRepEat. *bioRxiv* (2026).

Mosquito artwork thanks to Petra Korlevic.