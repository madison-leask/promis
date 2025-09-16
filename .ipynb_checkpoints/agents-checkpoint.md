Packaging PROMIS Pipeline for Bioconda
Overview of PROMIS and Packaging Objectives
PROMIS is a Snakemake-based bioinformatics pipeline designed for microsatellite instability (MSI) analysis using high-throughput sequencing data. It consists of a Snakemake workflow, several Python scripts, and reference CSV files (e.g. lists of MSI loci) that are used during analysis. Packaging PROMIS as a Conda package on Bioconda will greatly simplify installation and ensure reproducible setups for end-users. The goal is to create a Bioconda recipe for PROMIS so that users can install it via a single conda install command, instead of manually managing dependencies.
Bioconda packaging will target Linux (the primary environment for bioinformatics pipelines) and potentially macOS. Bioconda does not support Windows builds, so Windows can be ignored for now[1]. If the pipeline’s dependencies are available on macOS, the recipe can include macOS; otherwise, initial support can be Linux-only (with an option to extend to macOS later). We will focus on a noarch (architecture-independent) Python package if possible, since PROMIS is mostly Python and Snakemake (making it easier to support multiple platforms). The plan is to guide the Codex agent through preparing the code and writing the Bioconda recipe, while leaving low-level implementation details (exact syntax, minor version pinning, etc.) for the agent to decide.
Preparing the Codebase for Conda Packaging
Before writing the conda recipe, the PROMIS codebase should be organized for installation:
•	Versioning and Release: Determine a version number for PROMIS (e.g. 1.0.0 or 0.1.0). Bioconda requires a version in the recipe, so having a tagged release on a public repository (GitHub/GitLab) is helpful. Prepare a source tarball or Git tag URL for this version.
•	License and Metadata: Ensure the pipeline has an explicit open-source license. The recipe will need a license field (e.g. MIT, BSD-3, etc.) and ideally include a LICENSE file. Also note a short summary or description of the tool (e.g. "Snakemake pipeline for microsatellite instability detection"), and a homepage URL (repository URL) for the about section of the recipe[2].
•	Structure for Installation: If PROMIS is not already a Python package, consider reorganizing it for easier installation. A common approach is to create a Python module (e.g. a promis/ directory with __init__.py) and include the Snakefile and other resources within it. This way, the pipeline can be installed via pip and accessed as a module. Alternatively, if you prefer not to restructure heavily, you can still package the existing files by copying them into the conda environment in a known location (e.g. $PREFIX/share/promis). Either approach can work for Bioconda. Using a Python packaging approach (with a setup.py or pyproject.toml) is beneficial because Conda can then use pip install to handle installation[3]. This also makes it easier to create a console script entry point.
•	Data Files: PROMIS uses reference data CSV files (for MSI loci). These must be included in the package. If using Python packaging, list these files as package data in setup.py or include them via MANIFEST.in so that they get installed in the module’s directory. If using a manual install (copying files in build.sh), ensure to copy the CSV files to an appropriate directory in the environment (e.g. $PREFIX/share/promis/ or within the module) and adjust the pipeline to look for them at runtime. The pipeline’s code should not assume the data files are in the current working directory unless you plan to copy them there; it’s better to reference them relative to the installed package path.
•	Command-Line Interface: Plan how users will run the pipeline after installation. It’s best to provide an easy entry point, such as a shell script or Python console script named promis. This wrapper would invoke Snakemake on the installed Snakefile. Many Snakemake-based pipelines implement a CLI for ease of use – for example, the MeSS pipeline provides a mess command (mess run -i samples.tsv) to launch its Snakemake workflow[4]. Under the hood, this is achieved by using a small Python CLI (often with Click or Argparse) that calls Snakemake’s API or subprocess. Using such a CLI greatly improves user experience by allowing them to run the pipeline with a single command[5]. For PROMIS, you can create a script (e.g. promis in $PREFIX/bin) that does something like: snakemake -s /path/to/installed/Snakefile --configfile ... plus any arguments the user provides. If you package via setup.py, you can specify an entry point in setup.py (entry_points={'console_scripts': ['promis=promis.cli:main']}) to auto-install the promis command. If not, you can manually add a wrapper script in the conda recipe’s build steps. The key point is to automate pipeline execution for the user.
•	Gather Dependencies: Review the pipeline’s scripts and Snakefile to list all dependencies:
•	Python libraries: List all Python imports used in your scripts (e.g. pandas, biopython, pysam, pyyaml, etc.). These will go into the conda recipe under requirements: run. For example, if the scripts use Pandas and Biopython, those should be runtime deps (PROMIS might use Pandas for data frames or Biopython for sequence handling – similar pipeline recipes include those dependencies[6]). Check versions if certain features are needed, otherwise depend on a recent minimum version.
•	Snakemake: The pipeline obviously needs Snakemake to run. use snakemake version 9.9.0 Choose a version compatible with your Snakefile syntax.
•	External tools: List any command-line tools that the pipeline calls. For MSI analysis, PROMIS might use tools like msisensor-pro, msisensor2, MANTIS, or general tools such as samtools, bcftools, bedtools, etc., depending on implementation. Make sure each of these is available on Bioconda:
o	msisensor-pro – available on Bioconda as msisensor-pro[9] (if PROMIS uses it for MSI scoring).
o	samtools – available (samtools package).
o	htslib – if low-level usage (though samtools includes htslib).
o	bedtools or others if used.
o	If PROMIS uses an R script or library, list the R packages needed (these would also have conda packages, often prefixed with r-).
•	Ensure availability: Verify that each dependency is present either in Bioconda or Conda-Forge (the two channels used by Bioconda). If a needed dependency is missing, you have two choices: modify the pipeline to avoid it, or create a new conda recipe for that dependency as part of your Bioconda contribution[10]. (For example, if there's a specialized Python library not yet packaged, you might have to add it to Conda-Forge or Bioconda first). In most cases, common bioinformatics tools will already be present on Bioconda.
By the end of this preparation, you should have: a clear list of runtime dependencies, a decided approach for installation (pip vs manual copy), a defined entry point strategy, and a source archive for the pipeline code. With these in hand, you can proceed to writing the Bioconda recipe.
Writing the Bioconda Recipe (meta.yaml)
The Bioconda recipe for PROMIS will consist of a meta.yaml file (and possibly a simple build.sh script, if needed). Below are the key sections to include and how to fill them:
•	Package Metadata (package and about sections):
In meta.yaml, start by specifying the package name and version. Use lowercase for the name (e.g. promis). For example:

 	package:
  name: "promis"
  version: "1.0.0"
 	Under about:, provide the metadata:
•	summary: a one-line description of PROMIS (e.g. "Snakemake pipeline for detecting microsatellite instability from NGS data").
•	home: the URL to the project homepage or code repository (e.g. GitHub link).
•	license: the license name (e.g. MIT, BSD-3-Clause, etc.).
•	license_file: (optional) if a LICENSE text is included in the source, you can point to it so that it gets included in the package.
This mirrors what other pipeline recipes do – for instance, the zAMP pipeline recipe lists its homepage and license in this section[2].
•	Source Section (source):
Provide a source from which the conda build will fetch the code. Common options:
•	A URL to a tarball of the release (for GitHub, you can use an archive link of a tag, e.g. https://github.com/username/promis/archive/v1.0.0.tar.gz). Include the sha256 hash of the file for integrity.
•	Alternatively, you can use a Git URL and tag, but tarball + sha256 is simpler and preferred for fixed releases.
Make sure the source archive contains all necessary files (Snakefile, scripts, CSVs, etc.). If certain large files are excluded in the source distribution, adjust accordingly or fetch them separately (though for ~80KB CSV this is not an issue).
•	Build Section (build):
In this section, specify build-related settings:
•	noarch: python – since this is a pure-Python/Snakemake pipeline with no compiled C/C++ extensions, mark it as noarch. This tells conda the package is architecture-independent and will be built once (on Linux) and usable on all platforms[1] (except Windows, which is excluded by Bioconda policy). Noarch Python packages require the build and runtime to be flexible across Python versions; we’ll handle that in requirements.
•	Build Script or Commands: You have two approaches here:
a.	Use Pip (preferred if set up): If you prepared a setup.py or equivalent, you can use the pip install mechanism. For example:

 	build:
  noarch: python
  script: "{{ PYTHON }} -m pip install . -vv"
 	This will use Python’s pip to install the package from the source into the conda environment[3]. It assumes the source has a setup script. This approach automatically puts scripts into bin/ if defined as entry_points. It’s clean and leverages your Python packaging configuration.
b.	Manual Installation (if not using pip): If setting up a Python package is not feasible, you can write a simple build.sh (and include it in the recipe) to copy files. For example, the build script could:
c.	Create a directory in $PREFIX/share/promis and copy the Snakefile, config, CSV, and scripts there.
d.	Copy or create a launcher script in $PREFIX/bin/promis that calls Snakemake with the Snakefile. (Ensure to chmod +x it to be executable.)
e.	Alternatively, install the Snakefile and scripts under $PREFIX/bin or $PREFIX/lib/python3.X/site-packages as needed.
Using a build.sh gives you full control. In meta.yaml, you would then specify script: bash build.sh. However, if pip installation is possible, that is simpler and less error-prone.
•	Build Dependencies: Under requirements: host, list what’s needed to build the package. This typically includes:
o	python (with a version, e.g. >=3.8) – required in host to run the build (pip or any installation scripts).
o	pip – if you use the pip install method.
o	If using build.sh to do copying, you might not need anything beyond basic shell (which is always available) and snakemake is not needed at build time, only run time.
o	If you were compiling code, compilers would go here, but for this pipeline, likely not needed.
•	No compilation is expected, so no bld.bat (Windows build script) is needed – we explicitly skip Windows anyway.
•	Requirements Section (requirements: run):
This is one of the most important parts: it lists all runtime dependencies that must be present when someone installs PROMIS. Based on the earlier dependency gathering:
•	Python: specify the Python version range. For example, if you want to support Python 3.9 and above, you could put python >=3.9. It’s common to exclude Python 2 (which is EOL) and possibly limit to Python 3.11 or 3.10 if you know your code needs a certain version.
•	Snakemake: include snakemake-minimal. You can specify a broad version or pin to a major version. For instance: snakemake-minimal >=7.0 (or snakemake-minimal 7.* if you only tested on v7, or 8.* for v8). Many Bioconda recipes pin Snakemake’s version range to avoid future incompatible changes[11][7]. For a new pipeline, using the latest Snakemake is fine, but you might still restrict to the latest major (e.g., >=7,<8.999 if you want all 7.x and 8.x).
•	Workflow Tools: add all command-line tools used by the pipeline:
o	If PROMIS calls msisensor-pro, include msisensor-pro (and you might specify >=1.2 or a range if a particular version is required). Likewise for any other MSI detection tool.
o	If it uses samtools for BAM handling, add samtools (e.g. >=1.12 as used in similar pipelines[12]).
o	Include any other tools (e.g. bcftools, bedtools, python-bedtools if using pybedtools, etc.).
o	If R scripts are involved, list the R base and packages (e.g. r-base >=4.1 and specific r-<pkg> names).
•	Python Libraries: add all needed Python packages:
o	e.g. pandas, biopython, pysam, numpy, pyyaml, etc., with versions if needed. For instance, AQUAMIS lists numpy >=1.21, pandas >=1.3.5, pyyaml >=6 as dependencies[13][14]. Include whatever PROMIS requires.
o	If the pipeline uses a CLI framework like Click or Argparse, include those. (For example, if you implement the CLI with Click, add click >=8.1.)
o	The snaketool-utils library (and related, like attrmap or metasnek) – these appear in some pipeline recipes (MeSS, zAMP) to support their CLI[7][15]. If you decide to use the Snaketool framework or borrowed code from it, include these packages accordingly. If not, you can ignore them. (They are not required unless you explicitly used them in PROMIS.)
•	Pinning and Channels: Bioconda automatically pulls packages from conda-forge (for general packages) and bioconda (for bio-specific). You generally do not specify channels in the recipe (the build system handles that)[16]. Just ensure the names and versions match what exists. Pin versions only if necessary to avoid conflicts. For example, it’s good practice to avoid overly strict pins unless the pipeline is known to break with other versions. Something like msisensor-pro >=1.2 or pandas >=1.3 is usually sufficient to ensure a minimum version[12]. If a dependency has breaking changes, you can upper-bound it.
•	By listing all these, when a user installs PROMIS, conda will pull all these packages (and their dependencies) automatically, providing a fully functional environment for the pipeline.
•	Test Section (test):
Bioconda recipes include a test to verify the installation was successful. This is not meant to be an exhaustive pipeline run (which could be too slow or require big data), but rather a simple check:
•	The easiest test for a pipeline with a CLI is to call the help/version. For example, you can have:

 	test:
  commands:
    - promis --help
 	This will run promis --help inside a fresh environment (with only the runtime deps installed) and check for an exit code 0. It ensures that the promis command is properly installed and can start. If your CLI returns 0 on --help and prints usage, this is a good smoke test.
•	Alternatively, if no CLI wrapper is provided, you could test an import: python -c "import promis" (if you have a Python module). But a CLI test is preferred for an end-user tool. For instance, other CLI tools test their command in recipes (if a tool has subcommands, sometimes they test toolname --version or --help). Ensure that running promis --help doesn’t actually execute a long workflow (it should ideally just print usage and exit).
•	If the pipeline has a small example built-in or a dry-run option, you could use that in tests. For example, promis --dry-run -c config.yaml if a config and dry-run mode exist, but this might require including a tiny test file. In most cases, a simple help/usage test is sufficient[17].
•	No external dependencies in tests: The test environment will not include things outside your run requirements unless explicitly listed under test: requires:. Don’t rely on network access or large files. Keep it self-contained. If needed, you can include a very small test file in the source and call a command on it, but if not trivial, skip that. The goal is just to ensure installation isn’t broken.
•	Note: If you need to do a more complex test that requires additional packages (not in run reqs) or data, you could use a run_test.sh script in the recipe[18]. However, for PROMIS, this is likely unnecessary.
•	Additional Recipe Fields:
You can include other fields as needed:
•	extra: can hold extra metadata (not usually needed).
•	If the build produces multiple outputs (subpackages), that can be defined too, but in this case PROMIS will be a single package.
•	Skip/Platform selectors: Since we only support Linux (and maybe macOS), you might add selectors to skip unsupported platforms. For example, to skip Windows explicitly (though Bioconda won’t try Windows builds at all), you could do:

 	build:
  skip: true  # [win] 
 	And if for some reason you wanted to skip macOS (if a critical dependency isn’t available on macOS), you could add skip: true # [osx]. Generally, if noarch: python is used, you should not need OS-specific skips for Linux/macOS (noarch will build on Linux and also cover macOS). Only use skip if you run into problems with macOS builds of a dependency. As an example, if msisensor-pro lacked a macOS build, you might restrict to linux by using a selector like # [linux] on that dependency or skipping osx entirely. But try to support both Linux and macOS if possible, since many Bioconda packages do (Bioconda can build noarch or separate builds for linux-64 and osx-64 as long as all deps exist on both).
•	Recipe Example Reference:
It may help to look at an existing Bioconda recipe for a Snakemake pipeline to guide formatting. For instance, the MeSS recipe (metagenomic simulator) or zAMP recipe can serve as templates. They show how the meta.yaml is structured and how dependencies are listed[7][15]. Notably, those recipes include entries like snakemake-minimal, various Python libs, and even snaketool-utils for CLI support. You can pattern PROMIS’s recipe similarly, adapting the dependencies for MSI-specific ones.
In summary, the meta.yaml will gather everything we planned: the source of the code, how to install it, what it needs to run, and how to test it. Once this file (and any auxiliary build/test scripts) is written, we are ready to build and submit the package.
Handling Platform Compatibility
As mentioned, Linux is the primary target for Bioconda packages. The PROMIS recipe should work on a 64-bit Linux build, which is the default. Bioconda’s build system will also attempt to build on macOS by default if the recipe is not marked noarch. If we use noarch: python, we only do one build (on Linux) and that package will be usable on macOS as well, provided that all dependencies are available for macOS.
•	macOS Support: Check that each dependency is available on macOS. Most pure Python libraries will be (they are noarch themselves or have Mac builds on conda-forge). Many bioinformatics tools have macOS builds on Bioconda or conda-forge. For example, samtools, bedtools, msisensor-pro (it appears msisensor-pro has builds for OSX ARM64 at least[19], and possibly x86_64), etc. If one of the needed tools is not available on macOS, you have a few options:
•	If the tool is critical and has no Mac build, you might restrict the package to Linux only (using a selector as discussed).
•	If the tool is optional in your pipeline, you could modify the workflow to skip that part on Mac or use an alternative.
•	Since the user explicitly mentioned possibly adding other OS support if not too much work, it’s worth trying to include macOS. Many pipelines do support macOS by virtue of being noarch Python or having dependencies built for Mac. So, aim for noarch: python which will make the package available on macOS automatically, unless we hit a snag.
•	You do not need to do anything special for different CPU architectures (x86_64 vs ARM) in a noarch package. Conda will handle choosing the correct dependency builds (e.g., on Apple M1, it will use osx-arm64 builds of dependencies if available). If you were not using noarch, Bioconda CI would build separate packages for each platform (Linux, OSX). But noarch simplifies that.
•	Windows: As noted, skip Windows entirely. Bioconda doesn’t build Windows packages[1], and most of the bioinformatics tools won’t run on native Windows anyway. If someone on Windows wants to use PROMIS, they would typically do so via WSL or a Docker container, but that is outside the scope of Bioconda. You do not need to mention Windows in the recipe (and if using noarch: python, Windows is automatically excluded since Bioconda channel is not configured for it).
•	Testing on Platforms: After packaging, it’s good to test the installation on both Linux and Mac (if available) to ensure everything works. Since this is a noarch package, a single build should work identically. But if, for example, a path issue arises on Mac due to case sensitivity or a dependency is missing, it will show up when running on Mac. Keep an eye out for those. The Bioconda CI will not explicitly run the tool on Mac (if noarch), but if you had a Mac-specific build, it would run the same tests on Mac.
In summary, by using a noarch Python package and only including cross-platform dependencies, we aim to support Linux (primary) and macOS with minimal extra effort. This aligns with Bioconda’s typical usage in HPC and development environments.
Submitting to Bioconda and Verification
Once the recipe is written, the next steps involve contributing it to Bioconda’s repository and ensuring it passes all checks:
1.	Fork and Branch: Fork the bioconda/bioconda-recipes GitHub repository if you haven’t already. Then create a new git branch for your recipe addition[20] (e.g., add-promis or similar). This isolates your changes.
2.	Add Recipe Files: Within your fork, add a new directory under recipes/ named promis (matching the package name). Place the meta.yaml inside this folder. If you use a separate build.sh or run_test.sh, place those in the same promis/ folder. The directory structure will look like:
 	bioconda-recipes/recipes/promis/meta.yaml
bioconda-recipes/recipes/promis/build.sh        (if needed)
bioconda-recipes/recipes/promis/run_test.sh    (if needed)
 	Typically only meta.yaml is required, plus build.sh if you opted for manual installation.
3.	Linting and Local Build (Optional but Recommended): Bioconda provides a linting tool to catch common issues. You can run bioconda-utils lint or simply rely on the automated lint. It’s also a good idea to do a local test build with conda-build before submission[21]. For example:
 	conda create -n test-build -c conda-forge -c bioconda conda-build
conda activate test-build
conda-build recipes/promis
 	This will attempt to build the recipe locally. If there are missing dependencies or mistakes in meta.yaml, you’ll get errors that you can fix before making the pull request. Local building is a quick way to iterate on the recipe (and you can inspect the built package to ensure files are in place, and even install it to test the promis command).
4.	Commit and Push: Once the recipe is ready and passes local tests, commit the changes in Git. Write a clear commit message like “Add PROMIS pipeline 1.0.0 (new recipe)”. Push the branch to your fork[22].
5.	Open a Pull Request: On GitHub, create a PR from your branch to the bioconda-recipes master branch[23]. The PR template will ask you to confirm some items (e.g., that you ran lint, that your tool isn’t already available elsewhere, etc.). Fill that in appropriately (PROMIS is a new tool, so that’s fine).
6.	Bioconda CI Build: After opening the PR, GitHub Actions (via Bioconda’s CI) will automatically run:
7.	Linting: It checks the meta.yaml format (common issues include missing license, too long summary, etc.).
8.	Test build: It will attempt to build the recipe in an isolated environment on Linux (and if not noarch, on Mac too).
9.	Mulled tests: Bioconda has a “mulled” container test that creates a minimal container with just your package’s runtime requirements and runs the test commands[17]. This ensures that your test: commands indeed work with only the listed deps. If any step fails, the CI will mark it and provide logs. Common fixes might be: adjusting a dependency name or version if it wasn’t found, adding a missing dependency, correcting a path in test, etc. You can push further commits to your branch to update the recipe and the CI will re-run. Don’t be discouraged by failures; iterative fixes are normal[24].
10.	Review and Merge: Bioconda requires that a member of the Bioconda core team review and merge your PR. You can ping the gitter channel or wait — usually someone will take a look once it passes CI. Since this is a new recipe, they will verify things like: proper naming, license inclusion, no forbidden content, etc. Make sure to address any feedback. Once it’s approved, it will get merged into the master branch of bioconda-recipes[25].
11.	Package Availability: After merge, the Bioconda build system will upload the package to the Anaconda cloud. Within ~30 minutes (sometimes faster), the package becomes installable by users from the bioconda channel[26]. You can announce to users that they can do conda install -c bioconda promis (and ensure they also have -c conda-forge in their config, as per Bioconda setup). It’s good to do a fresh install test in a clean environment to double-check everything is okay.
12.	Maintenance: Going forward, remember to update the recipe for new versions of PROMIS. Each new release will require bumping the version in meta.yaml and perhaps adjusting dependencies if they change. The process for updates is similar (open PR with version bump). If any of PROMIS’s dependencies get updated or deprecated, you may need to update the recipe accordingly over time.
Throughout this process, you leave the fine-grained decisions (like exact version constraints, minor code tweaks for installation) to the Codex agent implementing the instructions. For example, Codex might decide on the appropriate pinning for snakemake-minimal or how exactly to formulate the wrapper script — these are “small decisions” that the agent can infer from context and examples. The key is that the high-level plan and requirements are clearly stated, so the agent knows what outcome to achieve.
Conclusion and Additional Tips
By following this structured plan, the Codex agent (or any developer) can package the PROMIS pipeline as a Bioconda package with minimal friction:
•	We started by outlining the context and goals: packaging a Snakemake workflow for easy distribution, focusing on Linux (and macOS if possible) and skipping Windows[1].
•	We emphasized preparing the pipeline for packaging, including organizing the code and identifying all requirements (ensuring nothing is missed, as missing dependencies are a common pitfall in packaging).
•	We provided guidance on writing the conda recipe (meta.yaml) with all necessary fields — from source, build instructions, to a thorough list of runtime dependencies and a simple test — referencing best practices from existing Bioconda recipes (e.g., inclusion of snakemake-minimal[8] and other libraries).
•	The importance of an entry point for usability was highlighted, with suggestions to implement a CLI wrapper (drawing on how other pipelines use a CLI to launch Snakemake[5]). This will make PROMIS user-friendly once installed.
•	Platform considerations were addressed so that Codex knows to mark the package appropriately (noarch, no Windows) and handle macOS compatibility if applicable.
•	Finally, we walked through the contribution process to Bioconda, from forking the recipes repo to PR, CI, and merge[23][25], so the agent understands the workflow and can even automate parts of it if needed.
Throughout the writing of the recipe, it’s beneficial for Codex to refer to Bioconda documentation and similar recipes for reference. Key references have been cited, which Codex can consult for specific syntax or examples (e.g., how dependencies are formatted or how tests are structured in recipes). This ensures the guidance is not just theoretical but grounded in real examples[7][15].
Once the recipe is built and merged, PROMIS will be available to users via Bioconda. They will be able to create an environment and install PROMIS with a single command, with all its dependencies (Python, Snakemake, and bioinformatics tools) pulled in automatically. This significantly lowers the barrier to using the pipeline and aligns with reproducibility standards in computational biology.
References:
•	Bioconda Contribution Workflow (steps to add a recipe)[27][23]
•	Bioconda Guidelines (Python package recipes, noarch, channel usage, platform support)[3][1][16]
•	Example Bioconda Recipes for Snakemake pipelines: MeSS (metagenomic simulator) and zAMP (amplicon pipeline) showing dependency lists and usage of snakemake-minimal[7][15]
•	Snaketool CLI concept used in pipelines like MeSS (for easier Snakemake execution)[5]. This inspired the approach to include a promis command for the pipeline.
•	Dependency management in recipes: Aquamis recipe listing bioinformatics tool dependencies[12] and typical Python libraries (pandas, numpy, etc.)[13], illustrating how to specify versions.
By adhering to these guidelines and using the cited examples as models, Codex should be able to generate a correct and effective Bioconda recipe for PROMIS, and handle the packaging process smoothly. The result will be a conda-installable PROMIS pipeline, streamlining its deployment for all users. [1][5]
________________________________________
[1] [3] [10] [16] [17] [18] Guidelines for bioconda recipes — Bioconda documentation
https://bioconda.github.io/contributor/guidelines.html
[2] [11] [15] Package Recipe 'zamp' — Bioconda documentation
https://bioconda.github.io/recipes/zamp/README.html
[4] GitHub - metagenlab/MeSS: Snakemake pipeline for simulating shotgun metagenomic samples
https://github.com/metagenlab/MeSS
[5] MeSS and assembly_finder: a toolkit for in silico metagenomic ...
https://pmc.ncbi.nlm.nih.gov/articles/PMC11755095/
[6] [7] Package Recipe 'mess' — Bioconda documentation
https://bioconda.github.io/recipes/mess/README.html
[8] [12] [13] [14] Package Recipe 'aquamis' — Bioconda documentation
https://bioconda.github.io/recipes/aquamis/README.html
[9] [19] Package Recipe 'msisensor-pro' — Bioconda documentation
https://bioconda.github.io/recipes/msisensor-pro/README.html
[20] [21] [22] [23] [24] [25] [26] [27] Contribution Workflow — Bioconda documentation
https://bioconda.github.io/contributor/workflow.html
