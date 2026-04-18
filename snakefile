import os
import sys
import yaml
from glob import glob


def strip_alignment_extension(path):
    filename = os.path.basename(path)
    for ext in (".bam", ".cram"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename

# Load default configurations from `config.yaml`
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

# Update `config` with command-line overrides
for arg in sys.argv[1:]:
    if "=" in arg:
        key, value = arg.split("=", 1)
        config[key] = value

# Extract required configuration values
repeats = config.get("repeats")
output_dir = config.get("output_dir")
scripts_dir = config.get("scripts_dir")
msi_deviation = float(config.get("msi_deviation", 0.1))
cytoband = config.get("cytoband")
min_reads = int(config.get("min_reads", 100))
bq_threshold = float(config.get("bq_threshold", 28))
mq_threshold = int(config.get("mq_threshold", 40))
min_dev_reads = int(config.get("min_dev_reads", 30))
call_by = config.get("call_by", "both")
min_dev_percent = float(config.get("min_dev_percent", 10.0))
min_length_percent = float(config.get("min_length_percent", 5.0))
use_GMM = bool(config.get("use_GMM", True))
balance_tolerance = float(config.get("balance_tolerance", 0.01))
min_total_reads = int(config.get("min_total_reads", 50))

# Optional inputs
bam_files_str = config.get("bam_files", "")
input_dir = config.get("input_dir", "")

# Validate and collect alignment paths
if bam_files_str:
    bam_files = [x.strip() for x in bam_files_str.split(",")]
elif input_dir:
    bam_files = sorted(
        glob(os.path.join(input_dir, "**", "*.bam"), recursive=True)
        + glob(os.path.join(input_dir, "**", "*.cram"), recursive=True)
    )
else:
    raise ValueError("Please provide either 'bam_files' or 'input_dir' in the config.")

if not bam_files:
    raise ValueError("No BAM or CRAM files found or specified.")

# Validate mandatory keys
if not repeats:
    raise ValueError("Please specify the 'repeats' file path in the config.")
if not output_dir:
    raise ValueError("Please specify the 'output_dir' path in the config.")
if not scripts_dir:
    raise ValueError("Please specify the 'scripts_dir' path in the config.")
if not bam_files:
    raise ValueError("Please provide alignment files using '--config bam_files=<file1.bam>,<file2.cram>,...'")
if not cytoband:
    raise ValueError("Please specify the 'cytoband' path in the config.")

# Derive sample names from BAM/CRAM file paths.
sample_names = [strip_alignment_extension(bam) for bam in bam_files]

if len(sample_names) != len(set(sample_names)):
    raise ValueError("Detected duplicate sample names after stripping .bam/.cram extensions.")

rule all:
    input:
        expand("{output_dir}/combined_results.csv", output_dir=output_dir, sample=sample_names),
        #expand("{output_dir}/{sample}/{sample}_barplot_MSI.pdf", output_dir=output_dir, sample=sample_names),
        #expand("{output_dir}/{sample}/{sample}_scatter_plot.pdf", output_dir=output_dir, sample=sample_names),
        #expand("{output_dir}/{sample}/{sample}_heatmap_plot.pdf", output_dir=output_dir, sample=sample_names),
        #expand("{output_dir}/{sample}/{sample}_cytoband_instability_plot.pdf", output_dir=output_dir, sample=sample_names),
        expand("{output_dir}/{sample}/{sample}_repeat_type_summary.csv", output_dir=output_dir, sample=sample_names),
        #expand("{output_dir}/{sample}/{sample}_repeat_unit_instability_barplot.pdf", output_dir=output_dir, sample=sample_names),
        #expand("{output_dir}/{sample}/{sample}_repeat_length_vs_instability_scatterplot.pdf", output_dir=output_dir, sample=sample_names)

rule extract_reads:
    input:
        bam=lambda wildcards: bam_files[sample_names.index(wildcards.sample)],
        repeats=repeats
    output:
        extracted_csv="{output_dir}/{sample}/{sample}_extracted_reads.csv"
    conda:
        "envs/promis_rules.yaml"
    shell:
        """
        python {scripts_dir}/extract_MSI_sequences.py \
            -b {input.bam} \
            -r {input.repeats} \
            -o {output.extracted_csv} \
            --bq_threshold {bq_threshold} \
            --mq_threshold {mq_threshold} \
            --min_reads {min_reads}
        """

rule analyze_repeats:
    input:
        extracted_csv="{output_dir}/{sample}/{sample}_extracted_reads.csv",
        repeats=repeats
    output:
        repeats_csv="{output_dir}/{sample}/{sample}_repeats_analysis.csv"
    conda:
        "envs/promis_rules.yaml"
    shell:
        """
        python {scripts_dir}/analyze_MSI_lengths.py -e "{input.extracted_csv}" -r "{input.repeats}" -o "{output.repeats_csv}" -t {msi_deviation}
        """

rule analyze_distribution:
    input:
        repeats_csv="{output_dir}/{sample}/{sample}_repeats_analysis.csv"
    output:
        distribution_csv="{output_dir}/{sample}/{sample}_distribution_analysis.csv"
    params:
        gmm_flag=lambda wildcards: "--use_GMM" if use_GMM else ""
    conda:
        "envs/promis_rules.yaml"
    shell:
        """
        python {scripts_dir}/analyze_MSI_distribution.py \
            -i {input.repeats_csv} \
            -o {output.distribution_csv} \
            --call_by {call_by} \
            --min_dev_reads {min_dev_reads} \
            --min_dev_percent {min_dev_percent} \
            --min_length_percent {min_length_percent} \
            {params.gmm_flag} \
            --balance_tolerance {balance_tolerance} \
            --min_total_reads {min_total_reads}
        """

rule plot_msi_status:
    input:
        distribution_csv="{output_dir}/{sample}/{sample}_distribution_analysis.csv"
    output:
        barplot="{output_dir}/{sample}/{sample}_barplot_MSI.pdf"
    conda:
        "envs/promis_rules.yaml"
    shell:
        """
        python {scripts_dir}/plot_MSI_results.py -i {input.distribution_csv} -o {output.barplot}
        """

rule plot_region_stats:
    input:
        distribution_csv="{output_dir}/{sample}/{sample}_distribution_analysis.csv"
    output:
        scatter_plot="{output_dir}/{sample}/{sample}_scatter_plot.pdf",
        heatmap_plot="{output_dir}/{sample}/{sample}_heatmap_plot.pdf",
        cytoband_plot="{output_dir}/{sample}/{sample}_cytoband_instability_plot.pdf"
    conda:
        "envs/promis_rules.yaml"
    shell:
        """
        python {scripts_dir}/plot_region_stats.py \
        -d {input.distribution_csv} \
        -s {output.scatter_plot} \
        -m {output.heatmap_plot} \
        -c {cytoband} \
        --cytoband_output {output.cytoband_plot}
        """
rule analyze_repeat_types:
    input:
        lengths_csv="{output_dir}/{sample}/{sample}_repeats_analysis.csv"
    output:
        summary_csv="{output_dir}/{sample}/{sample}_repeat_type_summary.csv",
        instability_plot="{output_dir}/{sample}/{sample}_repeat_unit_instability_barplot.pdf",
        length_vs_instability_plot="{output_dir}/{sample}/{sample}_repeat_length_vs_instability_scatterplot.pdf"
    conda:
        "envs/promis_rules.yaml"
    shell:
        """
        python {scripts_dir}/repeat_type_stats.py \
            -l {input.lengths_csv} \
            -o $(dirname {output.summary_csv}) \
            --plot_instability \
            --plot_length_vs_instability \
            --sample {wildcards.sample}
        """

rule combine_results:
    input:
        expand(
            "{output_dir}/{sample}/{sample}_distribution_analysis.csv",
            output_dir=output_dir,
            sample=sample_names,
        )
    output:
        combined_csv = "{output_dir}/combined_results.csv"
    run:
        import pandas as pd
        import re
        
        inputs = list(input)
        if not inputs:
            pd.DataFrame(columns=["Score", "Unstable regions", "Sample"]).to_csv(output.combined_csv, index=False)
            return
        
        combined_rows = []
        
        for result_file in inputs:
            df = pd.read_csv(result_file)
            summary_row = df[df["Chromosome"] == "Summary"]
            sample_name = result_file.split("/")[-2]
        
            if not summary_row.empty:
                status_str = str(summary_row.iloc[0]["MSI_Status"])
                print(f"Parsing: {status_str} from sample: {sample_name}")
                match = re.match(r"([\d.]+)%.*\((\d+)/(\d+)", status_str)
                if match:
                    score = float(match.group(1))
                    unstable_regions = f"{match.group(2)}/{match.group(3)}"
                    total_regions = int(match.group(3))
                    combined_rows.append({
                        "Score": score,
                        "Unstable regions": unstable_regions,
                        "Total regions": total_regions,
                        "Sample": sample_name
                    })
        
        combined_df = pd.DataFrame(combined_rows)
        combined_df.to_csv(output.combined_csv, index=False)
