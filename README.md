# PROMIS

PROMIS is a Snakemake-based workflow for detecting microsatellite instability
(MSI) events in next-generation sequencing data. The pipeline bundles the
required reference tables and analysis scripts so it can be launched with a
single command-line entry point.

## Installation

PROMIS is packaged for Conda via Bioconda. After configuring the Bioconda
channels, install it with:

```bash
conda install -c conda-forge -c bioconda promis
```

## Quick start

1. **Create a working configuration.** Use the packaged template as a starting
   point:

   ```bash
   promis --print-config > config.yaml
   ```

2. **Edit `config.yaml`.** Update the `output_dir` and provide either a
   comma-separated list of BAM files (`bam_files`) or an input directory
   (`input_dir`). The paths to the bundled repeat annotations and cytoband data
   are filled in automatically.

3. **Run the workflow.** Invoke Snakemake through the `promis` wrapper:

   ```bash
   promis -j 8 --configfile config.yaml --use-conda --keep-going
   ```

   Additional Snakemake options can be appended to the command (for example,
   `--config min_reads=10`). The `--use-conda` flag enables the provided rule
   environment specification (`femwell_msi.yaml`).

## Command-line interface

```
promis --help
```

Key options:

- `--print-config` – print the default configuration file and exit.
- `--workflow-dir` – show the location of the installed workflow files.
- `-j/--cores` – control parallelism.
- `--configfile` – run with a custom configuration YAML file.
- `--use-conda` – delegate rule execution to per-rule Conda environments.

## Packaged resources

The Conda package installs the Snakemake workflow, Python scripts, and reference
files under the `promis.workflow` Python module. They can be located
programmatically:

```python
import promis
print(promis.get_default_config_path())
print(promis.get_snakefile_path())
```

## Platform support

PROMIS is distributed as a `noarch` Python package targeting Linux, with macOS
support provided automatically when all dependencies are available on the
Bioconda and conda-forge channels.

## License

PROMIS is released under the [MIT License](LICENSE).
