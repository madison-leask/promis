import pandas as pd
import numpy as np
import argparse
import logging
import pickle

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

def load_model(model_path):
    """Load a trained GMM model from disk."""
    with open(model_path, "rb") as f:
        return pickle.load(f)

def apply_gmm_to_sample(sample_file, model, output_file):
    """Apply a trained GMM model to classify MSI status for a sample."""
    
    df = pd.read_csv(sample_file)

    # Filter only reads where Context_Match == "Pass"
    df = df[df["Context_Match"] == "Pass"]

    # Ensure required columns exist
    if "Total_Length_With_Extensions" not in df.columns or "Expected_Length" not in df.columns:
        logger.warning(f"Skipping {sample_file} - missing required columns!")
        return

    # Convert 'Total_Length_With_Extensions' to numeric, dropping errors
    df["Total_Length_With_Extensions"] = pd.to_numeric(df["Total_Length_With_Extensions"], errors="coerce")
    
    # Drop NaN values caused by invalid entries
    df = df.dropna(subset=["Total_Length_With_Extensions"])

    results = []
    grouped = df.groupby(["Chromosome", "Region_Start", "Region_End"])
    
    for (chrom, start, end), group in grouped:
        if (chrom, start, end) not in model:
            logger.warning(f"No trained GMM for {chrom}:{start}-{end}. Skipping.")
            continue

        total_lengths = group["Total_Length_With_Extensions"].values.reshape(-1, 1)
        gmm = model[(chrom, start, end)]
        gmm_means = gmm.means_.flatten()
        gmm_weights = gmm.weights_.flatten()

        # Identify dominant repeat length
        sorted_indices = np.argsort(gmm_weights)[::-1]
        dominant_length = round(gmm_means[sorted_indices[0]])

        # Compute deviation threshold (5% of dominant length)
        deviation_threshold = 0.05 * dominant_length

        # Identify deviating reads
        deviating_reads = total_lengths[
            (total_lengths < dominant_length - deviation_threshold) |
            (total_lengths > dominant_length + deviation_threshold)
        ]

        percent_deviating = (len(deviating_reads) / len(total_lengths)) * 100
        msi_status = "Unstable" if len(deviating_reads) > 5 else "Stable"

        # Compute additional statistics
        mean_length = np.mean(total_lengths)
        median_length = np.median(total_lengths)
        std_length = np.std(total_lengths)
        min_length = np.min(total_lengths)
        max_length = np.max(total_lengths)
        expected_length = group["Expected_Length"].iloc[0]

        results.append({
            "Chromosome": chrom,
            "Region_Start": start,
            "Region_End": end,
            "Expected_Length": expected_length,
            "GMM_Length": dominant_length,
            "Mean": mean_length,
            "Median": median_length,
            "StdDev": std_length,
            "Min": min_length,
            "Max": max_length,
            "% Deviating Reads": percent_deviating,
            "Deviating_Reads": len(deviating_reads),
            "MSI_Status": msi_status,
            "Expected_Length_Reads": len(total_lengths)
        })

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Calculate overall summary stats
    total_regions = len(results_df)
    unstable_regions = (results_df["MSI_Status"] == "Unstable").sum()
    percent_unstable = (unstable_regions / total_regions) * 100 if total_regions > 0 else 0

    # Add summary row
    summary_row = {
        "Chromosome": "Summary",
        "Region_Start": None,
        "Region_End": None,
        "Expected_Length": None,
        "GMM_Length": None,
        "Mean": None,
        "Median": None,
        "StdDev": None,
        "Min": None,
        "Max": None,
        "% Deviating Reads": None,
        "Deviating_Reads": None,
        "MSI_Status": f"{percent_unstable:.2f}% Unstable ({unstable_regions}/{total_regions} regions)",
        "Expected_Length_Reads": None
    }

    results_df = pd.concat([results_df, pd.DataFrame([summary_row])], ignore_index=True)

    # Ensure Chromosome is sorted naturally: chr1, chr2, ..., chrX, chrY
    def chromosome_sort_key(chrom):
        try:
            return int(chrom.replace("chr", ""))
        except ValueError:
            return float('inf') if chrom == "chrX" else float('inf') + 1

    results_df["Chromosome"] = results_df["Chromosome"].astype(str)  # Ensure string format
    results_df = results_df.sort_values(by=["Chromosome"], key=lambda x: x.map(chromosome_sort_key))

    # Save results
    results_df.to_csv(output_file, index=False)
    logger.info(f"GMM-based MSI classification saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply a global GMM model for MSI classification.")
    parser.add_argument("--input", required=True, help="Path to repeat length analysis CSV file.")
    parser.add_argument("--model", required=True, help="Path to the trained GMM model.")
    parser.add_argument("--output", required=True, help="Path to save the MSI classification results.")

    args = parser.parse_args()

    global_model = load_model(args.model)
    apply_gmm_to_sample(args.input, global_model, args.output)
