##  arg-rank-proj

## Overview

This repository contains code used to run MetaCompare2.0 on assembled metagenomic contigs and to generate all figures included in the manuscript. The workflow is fully reproducible and uses Snakemake, conda, and R.

# 1. Requirements

## Software

- MetaCompare2.0 (external tool)  
  Source: https://github.com/mrumi/MetaCompare2.0/  
  Version used: `[insert version]`

## Conda Environment

The conda environment contains:

- Diamond
- Prodigal
- mmseqs2
- Python
- Biopython
- pandas

## R (≥ 4.0) Packages

Required R packages:

- tidyverse
- vegan
- cowplot
- patchwork
- ggsignif
- ggpubr
- broom
- lme4
- readxl
- ggrain
- scales

# 2. Input Structure

The workflow expects the following directory layout:

```text
rawdata/
    sample1.fna
    sample2.fna
    ...

cip3_samples
```

Where:

- `rawdata/` contains assembled metagenomic contigs in FASTA format (`.fna`)
- `cip3_samples` is a plain text file containing one sample name per line

Example:

```text
sample1
sample2
sample3
```

# 3. Running MetaCompare with Snakemake

The `Snakefile` automates MetaCompare across all samples listed in `cip3_samples`.

## Run the Workflow

```bash
snakemake -s metacompare.smk \
    --cores 48 \
    --jobs 12 \
    --use-conda \
    -p -k
```

## Workflow Actions

This command will:

- Create per-sample output directories
- Run MetaCompare on each sample
- Log output to `logs/<sample>_mc.log`
- Produce `<sample>_out.txt` in `output/<sample>/`

# 4. Assembly Statistics

Assembly statistics were generated using the provided Python utilities:

```text
assembly_stats.py
combine_assembly_stats.py
```

## Example Usage

```bash
python assembly_stats.py rawdata/*.fna --out assembly_stats.tsv

python combine_assembly_stats.py assembly_stats/*.tsv \
    --out combined_stats.tsv
```

# 5. Figure Generation

All manuscript figures (`METFIG1–METFIG5`) are produced using:

```text
MCP_FIGURES.R
```

## Run

```bash
Rscript MCP_FIGURES.R
```

## Outputs

```text
METFIG1.png
METFIG2.png
METFIG3.png
METFIG4.png
METFIG5.png
```

# 6. Reproducibility Notes

- All code is provided exactly as used in the analysis.
- MetaCompare2.0 itself is not included in this repository; users should install it from the official repository:

  https://github.com/mrumi/MetaCompare2.0/

- Intermediate processed files (e.g., ARG count matrices) are included where necessary for reproducing figures.
- No cluster-specific commands or manual steps are required beyond running Snakemake and the R script.
