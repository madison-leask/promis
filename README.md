# INTRA-MSI Workflow

This repository contains a Snakemake pipeline for assessing microsatellite instability (MSI) in sequencing data.

## Usage

### 1. Copy the repository and install Snakemake
Clone this repository and ensure Snakemake is available. A convenient option is to
install it via Conda:

```bash
git clone https://github.com/your_org/INTRA-MSI.git
cd INTRA-MSI
conda install -c bioconda -c conda-forge snakemake
```

The pipeline automatically creates its own Conda environments when invoked with
`--use-conda`, so no manual environment creation is required.

### 2. Configure paths
Edit `config.yaml` so that all paths reflect your local setup. Important keys include:

- `input_dir` – directory containing input BAM files
- `output_dir` – where results will be written
- `scripts_dir` – path to this repository's scripts directory
- `repeats` – location of the MSI loci table
- `cytoband` – path to the cytoband annotation file
- `min_reads` – minimum deduplicated reads per locus (default: 100)
- `min_dev_reads` – minimum deviating reads to call instability (default: 30)
- `bq_threshold` – base quality cutoff for read extraction (default: 28)
- `mq_threshold` – mapping quality cutoff for read extraction (default: 40)
- `call_by` – how to call instability (`count`, `percent`, or `both`)
- `min_dev_percent` – percentage of deviating reads to call instability (default: 10.0)
- `min_length_percent` – minimum reference coverage percentage (default: 5.0)
- `use_GMM` – apply Gaussian Mixture Model (default: true)
- `balance_tolerance` – allowed imbalance for GMM components (default: 0.01)

### Logging

Analysis scripts such as `scripts/analyze_MSI_distribution.py` support
additional logging controls:

- `--info` – enable info-level logging
- `-v/--verbose` – enable debug-level logging

### 3. Run the workflow
After adjusting the configuration, execute Snakemake. The example below shows a
typical command, where `-np` performs a dry run. Omit `-np` to run the
pipeline:

```bash
snakemake --config \
  input_dir="/home/isilon/humangenetik/ANALYSES/Georgios/TCGA_MSI/data/1st_run_COAD/bam/" \
  output_dir="output/TCGA/COAD_250620" \
  -c all -j 16 --use-conda --keep-going --rerun-incomplete \
  -s snakefile_test -np
```

## Expected outputs
Results are generated inside `output_dir` as defined in the configuration. For each sample the workflow produces:

- `<sample>_extracted_reads.csv`
- `<sample>_repeats_analysis.csv`
- `<sample>_distribution_analysis.csv`
- PDF visualisations:
  - `<sample>_barplot_MSI.pdf`
  - `<sample>_scatter_plot.pdf`
  - `<sample>_heatmap_plot.pdf`
  - `<sample>_cytoband_instability_plot.pdf`
  - `<sample>_repeat_unit_instability_barplot.pdf`
  - `<sample>_repeat_length_vs_instability_scatterplot.pdf`
- A combined summary file: `combined_results.csv`
- 