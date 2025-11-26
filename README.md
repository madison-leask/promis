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
```bash
# Clone the repository
 git clone https://github.com/promis-bio/promis.git
 cd promis

# Create and activate the PROMIS environment
 conda env create -f promis/workflow/environment.yml
 conda activate promis
```

## Inputs
- Coordinate-sorted, indexed BAM files for each tumor sample
- Reference genome FASTA with accompanying index files (FAI, BWA/Bowtie2 if applicable)
- MSI loci metadata provided in `promis/workflow/database/`
- Optional sample sheet (columns: `sample`, `bam`) if you prefer to manage inputs outside the config file

## Configuration
1. Copy the template configuration:
   ```bash
   cp promis/workflow/config.yaml config.yaml
   ```
2. Edit `config.yaml` to set:
   - `output_dir`: destination for per-sample folders and combined results
   - `bam_files` **or** `input_dir`: explicit comma-separated BAM list or a directory to search recursively
   - `repeats`, `cytoband`, `scripts_dir`: override only if using custom resources
   - Thresholds such as `min_reads`, `min_dev_reads`, `bq_threshold`, `mq_threshold`, `min_dev_percent`, and `use_GMM`

Comments in the template describe each option. All bundled resource paths resolve automatically when left at their defaults.

## Running the pipeline
From the repository root (or any working directory containing your edited `config.yaml`), run:
```bash
snakemake --use-conda --cores 8 --configfile config.yaml
```
This executes the full PROMIS workflow, creating rule-specific conda environments from `promis/workflow/environment.yml` and writing outputs under `output_dir`.

## Outputs
Results are organized under `output_dir` (default `results/promis`):
- `<sample>/<sample>_extracted_reads.csv`: filtered reads spanning each MSI locus
- `<sample>/<sample>_repeats_analysis.csv`: inferred repeat lengths per locus
- `<sample>/<sample>_distribution_analysis.csv`: stability calls and MSI status per locus
- `<sample>/<sample>_barplot_MSI.pdf`: MSI status barplot
- `<sample>/<sample>_scatter_plot.pdf`, `<sample>/<sample>_heatmap_plot.pdf`, `<sample>/<sample>_cytoband_instability_plot.pdf`: region-level visualizations
- `<sample>/<sample>_repeat_type_summary.csv` plus accompanying instability plots
- `combined_results.csv`: cohort-level table summarizing MSI scores and unstable region counts

Key columns in `combined_results.csv`:
- `Sample`: sample identifier (from BAM basename)
- `Score`: MSI score (% unstable loci)
- `Unstable regions`: number of unstable loci / total loci assessed
- `Total regions`: total loci assessed

## Citation
Vlachos et al., PROMIS: tumor-only profiling of microsatellite instability, bioRxiv 2025, DOI: TBD

## Repository summary
- Updated workflow orchestration (`promis/workflow/Snakefile`) with clearer configuration handling and rule descriptions
- Consolidated conda environment definition at `promis/workflow/environment.yml`
- Refined default configuration (`promis/workflow/config.yaml`) and packaging metadata
- Removed pipeline cache/checkpoint files and obsolete environment stubs

Top-level structure:
- `Snakefile`: entry point that includes the packaged workflow
- `promis/`: Python package with workflow, scripts, and bundled reference data
- `promis/workflow/`: Snakefile, config template, environment file, and MSI reference tables
- `docs/`: design notes
- `conda/`: bioconda recipe scaffold

Known limitations (v1.0):
- Requires coordinate-sorted, indexed BAM inputs; CRAM is not yet supported
- MSI loci provided for hg38; using other genomes requires supplying compatible loci/cytoband files
- External aligner/index availability is assumed; PROMIS does not perform alignment
