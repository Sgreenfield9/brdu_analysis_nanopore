import os
import numpy as np
import pandas as pd
import sys


# Used to get to the root directory
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
# Ensure repo root is on sys.path so plotting module imports work
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

# Import used to use combined dataset for BrdU % and smoothing
from plotting.M_phase_chromosome_plotting import prepare_dataframe, smooth_counts

def get_export_dir() -> str:
    """
    Determine the output directory for saving BrdU%'s.
     The output directory is set in our .env.
    """
    output_dir = os.getenv("BRDU_PCT_OUTPUT_M_THRESHOLD")
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "output", "BrdU_pct")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def compute_brdU_pct(df: pd.DataFrame, min_distance: int) -> pd.DataFrame:
    """
    This function mostly does the same thing as the plotting code.
    We will keep track of the chromosome and starting point.
    We will then calculate the BrdU% using the BrdU count/Nmod
    * 100. We are excluding all percentages that are 0. We will 
    clean up any inf/nan that the code will get. 
    """
    out = df.copy()

    # Matched plotting logic: per chromosome sort + drop the duplicates
    out = out.sort_values(["chrom", "start"]).drop_duplicates(subset=["chrom", "start"])

    # We get the percntage for BrdU_count / Nmod to get percentage
    # We exclude percentages that will equal zero. 
    out["BrdU_pct"] = np.where(
        out["Nmod"] > 0,
        100.0 * out["BrdU_count"] / out["Nmod"],
        np.nan
    )

    out.loc[out["BrdU_pct"] <= 0, "BrdU_pct"] = np.nan

    # Match plotting: smooth per chromosome (window=1000) on BrdU_pct
    out["BrdU_pct"] = (
        out.groupby("chrom", sort=False)["BrdU_pct"]
        .transform(lambda s: smooth_counts(s, window=1000))
    )

    # Optional thinning: keep at most one row per min distance per chromosome
    if min_distance > 0:
        keep_idx = []
        for chrom, chrom_df in out.groupby("chrom", sort=False):
            last_kept = None
            for idx, start in zip(chrom_df.index, chrom_df["start"]):
                if last_kept is None or (start - last_kept) >= min_distance:
                    keep_idx.append(idx)
                    last_kept = start
        out = out.loc[keep_idx]

    # Clean up infinites/nan
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["chrom", "start", "end", "BrdU_pct"])

    return out

def export_sorted_brdU_pct_csv(
    df_with_pct: pd.DataFrame,
    export_dir: str,
    filename: str = "BrdU_pct_sorted_desc.csv",
    threshold_pct: float = 75.0,
):
    """
    Sort BrdU% highest -> lowest and export to CSV with chromosome + base range.
    """
    out = df_with_pct.copy()
    out = out[out["BrdU_pct"] <= threshold_pct]
    out = out.sort_values("BrdU_pct", ascending=False)

    # The columns we are keeping after parsing
    cols = ["chrom", "start", "end", "BrdU_pct", "BrdU_count", "Nmod"]
    cols = [c for c in cols if c in out.columns]
    out = out[cols]

    save_path = os.path.join(export_dir, filename)
    out.to_csv(save_path, index=False)
    print(f"Saved CSV: {save_path}")


def main():
    export_dir = get_export_dir()

    # Combined positive + negative dataframe
    df = prepare_dataframe()

    for min_distance in (2500, 5000, 10000):
        df_pct = compute_brdU_pct(df, min_distance=min_distance)
        suffix_kb = min_distance / 1000.0
        filename = f"BrdU_pct_sorted_desc_{suffix_kb:g}kb.csv"
        export_sorted_brdU_pct_csv(df_pct, export_dir, filename=filename)


if __name__ == "__main__":
    main()
