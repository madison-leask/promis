"""
===============================================================================
extract_MSI_sequences.py
===============================================================================

Description:
    This script extracts high-quality reads from a BAM file that fully span
    predefined microsatellite (MSI) repeat regions. It supports consensus BAMs
    with UMIs and enables filtering by base quality, mapping quality, and 
    sequence ambiguity.

    If Unique Molecular Identifiers (UMIs) are present in the RX, UR, or UB 
    SAM tags, the script will:
      - Extract the UMI per read
      - Deduplicate reads at each MSI locus using (Chromosome, Repeat_Coordinates, UMI)
      - Retain the highest-confidence read based on:
            1. Mapping quality
            2. Mean base quality
            3. Read length

Features:
    • Supports configurable base quality and mapping quality thresholds
    • Optional filtering of reads containing ambiguous 'N' bases
    • Automatically resolves chromosome naming mismatches (e.g., 'chr1' vs '1')
    • Logs UMI tag usage and deduplication status
    • Outputs a deduplicated CSV table with MSI region read summaries

Intended Use:
    Designed for pre-processing consensus BAM files in MSI detection pipelines 
    that require accurate, per-molecule analysis of microsatellite repeat lengths.

===============================================================================
"""
import pysam
import pandas as pd
import logging
import os
import argparse
import statistics
from rich.progress import track

# Configure logging
logger = logging.getLogger(__name__)

def load_repeat_coordinates(repeats_file):
    logger.info(f"Loading repeat coordinates from: {repeats_file}")
    repeats_df = pd.read_csv(repeats_file)
    logger.info(f"Columns in repeats_df: {repeats_df.columns}")
    return repeats_df

def check_chr_format(bam_file, repeats_df):
    repeats_df["Chromosome"] = repeats_df["Chromosome"].astype(str)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_chromosomes = bam.references
    bam_uses_chr = bam_chromosomes[0].startswith("chr")
    repeats_uses_chr = repeats_df["Chromosome"].iloc[0].startswith("chr")
    if bam_uses_chr and not repeats_uses_chr:
        repeats_df["Chromosome"] = repeats_df["Chromosome"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    elif not bam_uses_chr and repeats_uses_chr:
        repeats_df["Chromosome"] = repeats_df["Chromosome"].str.replace("^chr", "", regex=True)
    return repeats_df

def calculate_base_quality_stats(qualities):
    if not qualities:
        return {"Mean_Quality": None, "Median_Quality": None}
    mean_quality = round(statistics.mean(qualities), 2)
    median_quality = round(statistics.median(qualities), 2)
    return {"Mean_Quality": mean_quality, "Median_Quality": median_quality}

def extract_reads_from_bam(bam_file, repeats_df, output_file, bq_threshold, mq_threshold, keep_n, min_reads):
    try:
        sample_name = os.path.basename(bam_file).replace(".bam", "")
        logger.info(f"Processing sample: {sample_name}")
        bam = pysam.AlignmentFile(bam_file, "rb")
        extracted_data = []

        for _, row in track(
            repeats_df.iterrows(),
            total=len(repeats_df),
            description="Processing regions",
        ):
            chrom, start, end = row["Chromosome"], int(row["Start"]), int(row["End"])
            expected_repeat = row["Expected_Repeat"]
            repeat_start, repeat_end = int(row["Repeat_Start"]), int(row["Repeat_End"])
            logger.info(f"Processing region: {chrom}:{start}-{end}")
            reads = bam.fetch(contig=chrom, start=start, stop=end)
            for read in reads:
                qualities = read.query_qualities if read.query_qualities else []
                quality_stats = calculate_base_quality_stats(qualities)
                read_start = read.reference_start
                read_end = read.reference_start + read.query_length
                umi = None
                for tag in ("RX", "UR", "UB"):
                    try:
                        umi = read.get_tag(tag)
                        break
                    except KeyError:
                        continue
                if read_start <= repeat_start and read_end >= repeat_end:
                    if (quality_stats["Mean_Quality"] is not None and quality_stats["Mean_Quality"] > bq_threshold) and \
                       read.mapping_quality >= mq_threshold and \
                       (keep_n or "N" not in read.query_sequence):
                        extracted_data.append({
                            "Chromosome": chrom,
                            "Region_Start": start,
                            "Region_End": end,
                            "Read_Start": read_start,
                            "Read_End": read_end,
                            "Read_Name": read.query_name,
                            "Read_Sequence": read.query_sequence,
                            "UMI": umi,
                            "Mapping_Quality": read.mapping_quality,
                            "Read_Length": len(read.query_sequence) if read.query_sequence else 0,
                            "Expected_Repeat": expected_repeat,
                            "Repeat_Coordinates": f"{chrom}:{repeat_start}-{repeat_end}",
                            **quality_stats
                        })

        extracted_df = pd.DataFrame(extracted_data)
        grouped = extracted_df.groupby(["Chromosome", "Repeat_Coordinates"])
        final_rows = []

        for (chrom, repeat_coord), group in track(
            grouped,
            total=grouped.ngroups,
            description="Deduplicating",
        ):
            if "UMI" in group.columns and group["UMI"].notnull().any():
                dedup = group.sort_values(by=["Mapping_Quality", "Mean_Quality", "Read_Length"], ascending=False)
                dedup = dedup.drop_duplicates(subset=["Chromosome", "Repeat_Coordinates", "UMI"], keep="first")
            else:
                dedup = group

            if len(dedup) >= min_reads:
                final_rows.append(dedup)
                continue

            all_reads = []
            reads = bam.fetch(contig=chrom, start=int(group["Region_Start"].iloc[0]), stop=int(group["Region_End"].iloc[0]))
            for read in reads:
                read_start = read.reference_start
                read_end = read.reference_start + read.query_length
                if read_start > int(group["Read_Start"].iloc[0]) or read_end < int(group["Read_End"].iloc[0]):
                    continue
                qualities = read.query_qualities if read.query_qualities else []
                q_stats = calculate_base_quality_stats(qualities)
                if not keep_n and read.query_sequence and "N" in read.query_sequence:
                    continue
                umi = None
                for tag in ("RX", "UR", "UB"):
                    try:
                        umi = read.get_tag(tag)
                        break
                    except KeyError:
                        continue
                all_reads.append({
                    "Chromosome": chrom,
                    "Region_Start": int(group["Region_Start"].iloc[0]),
                    "Region_End": int(group["Region_End"].iloc[0]),
                    "Read_Start": read_start,
                    "Read_End": read_end,
                    "Read_Name": read.query_name,
                    "Read_Sequence": read.query_sequence,
                    "UMI": umi,
                    "Mapping_Quality": read.mapping_quality,
                    "Read_Length": len(read.query_sequence) if read.query_sequence else 0,
                    "Expected_Repeat": group["Expected_Repeat"].iloc[0],
                    "Repeat_Coordinates": repeat_coord,
                    **q_stats
                })

            fallback_df = pd.DataFrame(all_reads)
            supplemental = pd.DataFrame()
            if not fallback_df.empty:
                # Score and sort fallback reads by quality
                fallback_df["Quality_Score"] = fallback_df["Mapping_Quality"] + fallback_df["Mean_Quality"]
                fallback_df = fallback_df.sort_values(by="Quality_Score", ascending=False)
                
                # Calculate how many more reads needed
                n_missing = max(0, min_reads - len(dedup))
                
                # Take top-scoring reads
                supplemental = fallback_df.head(n_missing)

            # Append deduplicated and supplemental reads for this locus
            combined = pd.concat([dedup, supplemental], ignore_index=True)
            final_rows.append(combined)
                
        final_df = pd.concat(final_rows, ignore_index=True)
        final_df.to_csv(output_file, index=False)
        bam.close()
        logger.info(f"Deduplicated reads saved to: {output_file}")

    except Exception as e:
        logger.error(f"An error occurred while processing the BAM file: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads from BAM file for specific repeat coordinates.")
    parser.add_argument("-b", "--bam", required=True, help="Path to the BAM file.")
    parser.add_argument("-r", "--repeats", required=True, help="Path to the CSV file with repeat coordinates.")
    parser.add_argument("-o", "--output", required=True, help="Path to save the extracted reads CSV file.")
    parser.add_argument("--bq_threshold", type=float, default=38, help="Base quality threshold (default: 38)")
    parser.add_argument("--mq_threshold", type=int, default=58, help="Mapping quality threshold (default: 58)")
    parser.add_argument("--keep_n", type=lambda x: (str(x).lower() == 'false'), default=False, help="Whether to keep reads containing 'N' bases (default: False)")
    parser.add_argument("--min_reads", type=int, default=10, help="Minimum number of deduplicated reads per region (default: 10)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug-level logging.")
    parser.add_argument("--info", action="store_true", help="Enable info-level logging.")

    args = parser.parse_args()
    log_level = (
        logging.DEBUG if args.verbose
        else logging.INFO if args.info
        else logging.WARNING
    )
    logging.basicConfig(level=log_level, format='%(levelname)s:%(message)s')
    logger.setLevel(log_level)
    repeat_coords = load_repeat_coordinates(args.repeats)
    repeat_coords = check_chr_format(args.bam, repeat_coords)
    extract_reads_from_bam(args.bam, repeat_coords, args.output, args.bq_threshold, args.mq_threshold, args.keep_n, args.min_reads)

