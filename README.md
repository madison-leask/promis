# PROMIS: PROfiling of Microsatellite InStability

## Overview

PROMIS is a tumor-only, reference-free microsatellite instability (MSI) caller built on Snakemake. It models repeat-length mixtures of microsatellite loci and reports a continuous MSI score plus a binary MSI-high/MSS classification for whole-exome, whole-genome, targeted panel, on tissue or cell-free DNA sequencing data.

## Features

* Tumor-only workflow (no matched normal required)
* Discrete mixture modeling at microsatellite loci
* Continuous MSI score and binary MSI-high/MSS classification
* Works on exome, genome, targeted panels, and cfDNA
* Bundled reference loci, plotting scripts, and a reproducible execution setup

---

## Quickstart (recommended: run from the repository)

### 1) Clone the repository

```bash
git clone https://github.com/vlachosg37/promis.git
cd promis
```

### 2) Create the driver environment (Snakemake runner)

This environment runs Snakemake and creates rule-specific environments via `--use-conda`.

```bash
conda env create -f envs/promis_env.yaml
conda activate promis_env
```

Recommended (reproducibility and fewer solver issues):

```bash
conda config --env --set channel_priority strict
```

### 3) Prepare a configuration file

Edit paths for your setup:

```bash
nano config/config.yaml
```

At minimum, set:

At minimum, set:
- `input_dir`: directory containing coordinate-sorted BAMs and their indexes (`.bam` and matching `.bai`).
- `output_dir`: directory where PROMIS will write results (recommended: do not end with `/`).
- Any other reference/resource paths required by your run. By default, PROMIS uses `database/MSI_loci_hg38_exonic.csv`, which provides broad coverage suitable for WGS and WES. For targeted panels, see [Optional: Discovering microsatellite loci](#optional-discovering-microsatellite-loci).


### 4) Dry run (sanity check)

```bash
snakemake -np -c 1 --use-conda
```

### 5) Run

```bash
snakemake -c all --use-conda
```

---

## Notes on reproducibility

### Strict channel priority (important)

PROMIS assumes `conda-forge` + `bioconda` and **strict** channel priority. If your installation uses `defaults`, environment solves may fail or become non-reproducible.

Environment-local recommended settings:

```bash
conda config --env --set channel_priority strict
conda config --env --remove channels defaults 2>/dev/null || true
conda config --env --add channels conda-forge
conda config --env --add channels bioconda
```

### Clearing cached rule environments

If you change environment YAMLs or channels, clear Snakemake’s cached envs:

```bash
rm -rf .snakemake/conda
```

---

## Inputs

* Coordinate-sorted, indexed BAM (`.bam` + `.bai`) per tumor sample
* MSI loci table shipped in `database/` (or a custom loci file)
* Any additional resources specified in your `config.yaml`

## Outputs

Results are written under `output_dir`. Typical outputs include:

* `<sample>/<sample>_extracted_reads.csv`
* `<sample>/<sample>_repeats_analysis.csv`
* `<sample>/<sample>_distribution_analysis.csv`
* Per-sample PDFs (plots)
* `combined_results.csv` (cohort-level summary)

---

## Optional: Discovering microsatellite loci

PROMIS includes a helper to scan a reference genome for microsatellite loci. This is usefull if you have a targeted panel and you want to use loci found there. In general you 

```bash
python scripts/preprocess/find_MS_sites.py \
  --reference hg38.fa \
  --bed panel_targets.bed \
  --output database/panel_msi_loci.csv
```

## Use a representative BAM (coverage-restricted scanning)

```bash
python scripts/preprocess/find_MS_sites.py \
  --reference hg38.fa \
  --bam sample.bam \
  --min-coverage 30 \
  --output database/panel_msi_loci.csv
```

Use the resulting CSV as a custom loci file via your configuration.

---

## Citation

Vlachos et al., PROMIS: tumor-only profiling of microsatellite instability, bioRxiv (2025), DOI: TBD
