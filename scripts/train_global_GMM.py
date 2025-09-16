import pandas as pd
import numpy as np
import argparse
import logging
import pickle
from sklearn.mixture import GaussianMixture

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

def load_repeat_lengths(files):
    """Load repeat length data from multiple samples into a single DataFrame, filtering only 'Pass' reads."""
    all_data = []
    for file in files:
        logger.info(f"Loading data from: {file}")
        df = pd.read_csv(file)

        # Ensure required columns exist
        if "Total_Length_With_Extensions" in df.columns and "Expected_Length" in df.columns and "Context_Match" in df.columns:
            # **Filter only rows where Context_Match == "Pass"**
            df = df[df["Context_Match"] == "Pass"]

            # Convert 'Total_Length_With_Extensions' to numeric and drop invalid entries
            df["Total_Length_With_Extensions"] = pd.to_numeric(df["Total_Length_With_Extensions"], errors="coerce")
            valid_data = df.dropna(subset=["Total_Length_With_Extensions"])  # Drop rows with NaN

            if valid_data.empty:
                logger.warning(f"Skipping {file} - No valid reads after filtering!")
            else:
                all_data.append(valid_data[["Chromosome", "Region_Start", "Region_End", "Total_Length_With_Extensions", "Expected_Length"]])
        else:
            logger.warning(f"Skipping {file} - missing required columns!")

    return pd.concat(all_data, ignore_index=True) if all_data else None

def train_gmm(df):
    """Train a Gaussian Mixture Model (GMM) per microsatellite region."""
    gmm_models = {}

    grouped = df.groupby(["Chromosome", "Region_Start", "Region_End"])
    for (chrom, start, end), group in grouped:
        total_lengths = group["Total_Length_With_Extensions"].dropna().values.reshape(-1, 1)

        if len(total_lengths) < 10:  # Not enough data to train GMM
            logger.warning(f"Skipping {chrom}:{start}-{end} (too few valid reads).")
            continue

        gmm = GaussianMixture(n_components=2, random_state=42)
        gmm.fit(total_lengths)
        gmm_models[(chrom, start, end)] = gmm

    return gmm_models

def save_model(model, output_path):
    """Save trained GMM model to disk."""
    with open(output_path, "wb") as f:
        pickle.dump(model, f)
    logger.info(f"Trained global GMM model saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a global GMM model for MSI classification.")
    parser.add_argument("--input", nargs="+", required=True, help="Paths to repeat length analysis CSV files.")
    parser.add_argument("--output", required=True, help="Path to save the trained GMM model.")

    args = parser.parse_args()

    repeat_data = load_repeat_lengths(args.input)
    if repeat_data is not None:
        global_model = train_gmm(repeat_data)
        save_model(global_model, args.output)
    else:
        logger.error("No valid data available for training after filtering! Exiting.")
