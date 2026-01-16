import os
import numpy as np
import pandas as pd
import sys


# Used to get to the root directory
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
# Ensure repo root is on sys.path so plotting module imports work
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

# Import used to use combined dataset for BrdU %
from plotting.M_phase_chromosome_plotting import prepare_dataframe

def get_export_dir() -> str:
    """
    Determine the output directory for saving BrdU%'s
    """
    output_dir = os.getenv("BRDU_PCT_OUTPUT_M")
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "output", "BrdU_pct")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def compute_brdU_pct(df: pd.DataFrame) -> pd.DataFrame:
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

    # Clean up infinites/nan
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["chrom", "start", "end", "BrdU_pct"])

    return out

def export_sorted_brdU_pct_csv(df_with_pct: pd.DataFrame, export_dir: str, filename: str = "BrdU_pct_sorted_desc.csv"):
    """
    Sort BrdU% highest -> lowest and export to CSV with chromosome + base range.
    """
    out = df_with_pct.sort_values("BrdU_pct", ascending=False).copy()

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

    df_pct = compute_brdU_pct(df)
    export_sorted_brdU_pct_csv(df_pct, export_dir)


if __name__ == "__main__":
    main()
