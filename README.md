# PROMIS: PROfiling of Microsatellite InStability

## Overview
PROMIS is a tumor-only, reference-free microsatellite instability (MSI) caller built on Snakemake. It models repeat-length mixtures at hundreds of microsatellite loci and reports continuous MSI scores plus binary MSI-high/MSS status for whole-exome, whole-genome, targeted panel, and cell-free DNA sequencing data.

## Features
- Tumor-only workflow (no matched normal required)
- Discrete mixture modeling at microsatellite loci
- Continuous MSI score and binary MSI-high/MSS classification
- Works on exome, genome, and targeted panel data (including cfDNA)
- Bundled reference loci, plotting scripts, and a reproducible conda environment

## Installation

### Via conda (recommended)
PROMIS is available on Bioconda.

```bash
conda install -c bioconda -c conda-forge promis
```

### From source

```bash
git clone https://github.com/promis-bio/promis.git
cd promis

conda env create -f promis/workflow/environment.yml
conda activate promis
```

## Running the pipeline
Use the packaged CLI wrapper (defaults shown):

```bash
# From a directory containing your config.yaml:
promis --cores 8

# Or, explicitly:
promis --configfile config.yaml --cores 8
```

By default PROMIS:
- looks for `config.yaml` in the directory where you run the command,
- executes the packaged workflow, and
- enables `--use-conda` unless disabled.

Additional Snakemake options (for example `--config input_dir=... output_dir=...`) can be passed through after the known PROMIS options.

## Configuration
Place your `config.yaml` in the directory where you invoke `promis`, or point `--configfile` to its location. A template is packaged at `promis/workflow/config.yaml`; copy and edit it to set:
- `output_dir`: destination for per-sample folders and combined results
- `bam_files` **or** `input_dir`: explicit comma-separated BAM list or a directory to search recursively
- `repeats`, `cytoband`, `scripts_dir`: override only if using custom resources
- Thresholds such as `min_reads`, `min_dev_reads`, `bq_threshold`, `mq_threshold`, `min_dev_percent`, and `use_GMM`

Comments in the template describe each option. All bundled resource paths resolve automatically when left at their defaults.

## Inputs
- Coordinate-sorted, indexed BAM files for each tumor sample
- Reference genome FASTA with accompanying index files (FAI, BWA/Bowtie2 if applicable)
- MSI loci metadata provided in `promis/workflow/database/`
- Optional sample sheet (columns: `sample`, `bam`) if you prefer to manage inputs outside the config file

## Outputs
Results are organized under `output_dir` (default `results/promis`):
- `<sample>/<sample>_extracted_reads.csv`: filtered reads spanning each MSI locus
- `<sample>/<sample>_repeats_analysis.csv`: inferred repeat lengths per locus
- `<sample>/<sample>_distribution_analysis.csv`: stability calls and MSI status per locus
- `<sample>/<sample>_barplot_MSI.pdf`: MSI status barplot
- `<sample>/<sample>_scatter_plot.pdf`, `<sample>/<sample>_heatmap_plot.pdf`, `<sample>/<sample>_cytoband_instability_plot.pdf`: region-level visualizations
- `<sample>/<sample>_repeat_type_summary.csv` plus accompanying instability plots
- `combined_results.csv`: cohort-level table summarizing MSI scores and unstable region counts

## Optional: Discovering microsatellite loci
PROMIS ships with a helper CLI, `promis-find-ms-sites`, to scan a reference genome for microsatellite loci:

```bash
promis-find-ms-sites \
  --reference hg38.fa \
  --output hg38_msi_loci.csv \
  --bam example.bam \
  --min-coverage 30
```

The resulting CSV can be used as a custom loci file in the PROMIS configuration.

## Citation
Vlachos et al., PROMIS: tumor-only profiling of microsatellite instability, bioRxiv (2025), DOI: TBD
