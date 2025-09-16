# Architecture Decision Record: Package PROMIS workflow for Bioconda

## Status
Accepted

## Context
PROMIS previously existed as a flat Snakemake repository that required manual
path management and direct invocation of the Snakefile. The goal is to publish a
Bioconda package, which demands predictable installation paths, explicit runtime
dependencies, and a user-friendly entry point. The repository also bundled
reference CSV files and Python scripts that needed to be accessible after
installation without assuming a cloned Git checkout.

## Decision
Package the workflow as a Python module that exposes a console entry point and
includes all workflow assets as installable data. Provide Python packaging
metadata (pyproject.toml) and a Bioconda recipe that installs the project via
`pip`. Move the Snakefile, configuration template, rule environment file,
reference tables, and helper scripts into the `promis.workflow` package so they
can be resolved through `importlib.resources`. Introduce a `promis` CLI wrapper
that shells out to Snakemake with sensible defaults and helper options such as
`--print-config` and `--workflow-dir`.

### Technical Approach
- Added a `promis` Python package with `__init__.py`, `workflow/__init__.py`, and
a CLI module (`promis/cli.py`).
- Refactored the Snakefile to locate packaged resources relative to its module
path and to reference the bundled Conda environment file through an absolute
variable (`ENV_FILE`).
- Updated the default configuration to remove hard-coded script paths and to
work with the packaged resource layout.
- Declared project metadata in `pyproject.toml`, defined package data patterns in
`MANIFEST.in`, and registered the `promis` console script.
- Authored a Bioconda recipe (`conda/meta.yaml`) that builds the package in
`noarch: python` mode by invoking `python -m pip install . -vv` and declares all
Python runtime dependencies (pandas, numpy, scikit-learn, matplotlib-base,
seaborn, pysam, pyyaml, rich, tqdm, numba, snakemake-minimal).
- Installed the workflow assets into the package and updated the README with
Bioconda installation and CLI usage instructions.

### Key Components
- `promis/cli.py` – command-line entry point that invokes Snakemake.
- `promis/workflow/Snakefile` – Snakemake workflow adjusted to packaged paths.
- `pyproject.toml` and `MANIFEST.in` – Python packaging metadata and data
inclusion rules.
- `conda/meta.yaml` – Bioconda recipe describing build strategy and dependencies.

## Consequences

### Positive
- Users can install PROMIS through `conda install -c bioconda promis` with all
runtime dependencies resolved automatically.
- The `promis` CLI provides a single command to run the workflow and discover
packaged resources, simplifying documentation and support.
- Packaging resources inside the module removes the need for manual path
configuration and allows Snakemake to access reference files regardless of the
working directory.

### Negative
- The package now mandates Python ≥3.9 to support `importlib.resources.files` and
the selected dependency versions.
- Large third-party dependencies (e.g., numba, pysam, scikit-learn) are required
at runtime, increasing installation size and potentially limiting platform
support if builds are unavailable.
- Contributors must update package metadata and the Bioconda recipe for every new
release.

### Risks
- macOS availability of heavy dependencies (e.g., numba, pysam) could delay
Bioconda publication or necessitate Linux-only selectors.
- The CLI depends on the external `snakemake` executable provided by
`snakemake-minimal`; mismatched versions could cause runtime failures.
- Bundled data must stay synchronised with the workflow; missing files would
break packaged installations.

## Alternatives Considered

### Alternative 1: Install resources into `$PREFIX/share/promis`
**Description**: Copy workflow files into a shared data directory during Conda
build and provide a shell wrapper that references them.
**Pros**: Avoids Python packaging overhead; mirrors many legacy Bioconda
pipelines.
**Cons**: Requires custom install logic in `build.sh`, lacks Python-level access
to resources, and complicates reuse from other Python tooling.

### Alternative 2: Keep repository un-packaged and document manual setup
**Description**: Continue distributing PROMIS as a Git repository that users
clone and run directly.
**Pros**: Minimal engineering effort and no packaging metadata to maintain.
**Cons**: Leaves dependency management to users, prevents Conda distribution,
and provides no standardized CLI entry point.

## Implementation Details

### Dependencies
- Python ≥3.9
- pandas, numpy, scikit-learn, matplotlib-base, seaborn, pysam, pyyaml, rich,
  tqdm, numba, snakemake-minimal

### Configuration
- Default configuration template stored at `promis/workflow/config.yaml`.
- Resource resolver `_resolve_data_path` in the Snakefile makes paths relative to
the packaged workflow directory.
- Rule environments use `promis/workflow/femwell_msi.yaml` with added rich/tqdm
packages.

### Performance Considerations
- No runtime optimisations introduced; performance characteristics are governed
by existing Snakemake rules and dependencies. Packaging ensures consistent
resource loading, which may reduce path-related errors but does not change
computational complexity.

## Validation
- `python -m pip install --no-deps -e .` to ensure editable installation works
with the new packaging metadata.
- `python -m promis.cli --help` and `python -m promis.cli --print-config` to
confirm the CLI runs and locates bundled assets.

## References
- `pyproject.toml`
- `conda/meta.yaml`
- `promis/cli.py`
- `promis/workflow/Snakefile`

## Review Schedule
Review this decision when preparing future PROMIS releases or when dependency
constraints change (e.g., if numba/snakemake updates require altering the Python
minimum version). Revisit if Bioconda build issues necessitate a different
installation layout.
