"""
Chromosome heatmap: rows = chromosomes, cols = genomic bins (kb),
color = BrdU% (weighted by coverage Nmod)

Input: same bedgraph dataset used for genome browser:
columns = chrom, start, end, frac_mod, Nmod
"""

import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dotenv import load_dotenv

# Mapping for genbank to chromosome
genbank_to_chr = {
    "CM007964.1": "1",  "CM007965.1": "2",  "CM007966.1": "3",  "CM007967.1": "4",
    "CM007968.1": "5",  "CM007969.1": "6",  "CM007970.1": "7",  "CM007971.1": "8",
    "CM007972.1": "9",  "CM007973.1": "10", "CM007974.1": "11", "CM007975.1": "12",
    "CM007976.1": "13", "CM007977.1": "14", "CM007978.1": "15", "CM007979.1": "16",
    "CM007980.1": "p2-micron", "CM007981.1": "MT"
}

# Length of our yeast chromosomes 
chrom_lengths = {
    "1": 230218, "2": 813184, "3": 316620, "4": 1531933,
    "5": 576874, "6": 270161, "7": 1090940, "8": 562643,
    "9": 439888, "10": 745751, "11": 666816, "12": 1078177,
    "13": 924431, "14": 784333, "15": 1091291, "16": 948066,
    "MT": 85779, "p2-micron": 6318
}


def map_chromosome(chrom_id: str) -> Optional[str]:
    if chrom_id in genbank_to_chr:
        return genbank_to_chr[chrom_id]
    return None


def chrom_sort_key(x: str):
    if x.isdigit():
        return (0, int(x))
    return (1, x)

# Load the bedgraph data (we used the same data for our genomic browsers)
REPO_ROOT = Path(__file__).resolve().parents[4]
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")


def load_bedgraph_data():
    pos_path = os.getenv("POSITIVE_BEDGRAPH_M")
    neg_path = os.getenv("NEGATIVE_BEDGRAPH_M")
    if not pos_path or not neg_path:
        raise ValueError("POSITIVE_BEDGRAPH_M and NEGATIVE_BEDGRAPH_M must be set in env/.env")

    pos_df = pd.read_csv(
        pos_path, 
        sep="\t", 
        header=None,
        names=["chrom", "start", "end", "frac_mod", "Nmod"])
    
    neg_df = pd.read_csv(
        neg_path, 
        sep="\t", 
        header=None,
        names=["chrom", "start", "end", "frac_mod", "Nmod"])
    return pos_df, neg_df


def prepare_dataframe():
    """
    Combining our negative and positive strand datasets. 
    Then, we will map the chromosome ID to the correct chromosome number
    We then, put all of our data into panda dataframes (easier to plot this way with matplotlib)
    """
    pos_df, neg_df = load_bedgraph_data()
    df = pd.concat([pos_df, neg_df], ignore_index=True)

    df["chrom"] = df["chrom"].map(map_chromosome)
    df = df.dropna(subset=["chrom"]).copy()
    df["chrom"] = df["chrom"].astype(str)

    # BrdU percent per interval (0–100)
    df["brdu_pct"] = 100.0 * df["frac_mod"].astype(float)

    # basic cleanup
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    df["Nmod"] = df["Nmod"].astype(float)

    # keep only positive coverage
    df = df[df["Nmod"] > 0].copy()
    return df


def get_output_dir():
    """
    Check if output directory exists in the env/.env file. If not create a new directory. 
    """
    out = os.getenv("OUTPUT_DIR_M_HEATMAP")
    if not out:
        out = os.path.join(os.getcwd(), "output", "M_phase_pileup", "heatmaps")
    os.makedirs(out, exist_ok=True)
    return out


# Heatmap builder
def make_chr_heatmap(
    df: pd.DataFrame,
    chroms: list[str],
    bin_size_bp: int = 1000,
    min_cov: float = 10.0,
    include_extra: bool = False
):
    """
    Returns:
      H: (n_chroms, n_bins_max) matrix of BrdU% (weighted mean by Nmod)
      x_kb: x-axis bin centers in kb (0..max)
      chroms_used: row labels
    """
    if not include_extra:
        chroms_used = [c for c in chroms if c.isdigit()]
    else:
        chroms_used = chroms[:]

    max_len = max(chrom_lengths[c] for c in chroms_used)
    n_bins = int(np.ceil(max_len / bin_size_bp))
    x_kb = (np.arange(n_bins) + 0.5) * (bin_size_bp / 1000.0)

    H = np.full((len(chroms_used), n_bins), np.nan, dtype=float)

    df = df.copy()
    df["bin"] = (df["start"] // bin_size_bp).astype(np.int64)

    for r, chrom in enumerate(chroms_used):
        chr_len = chrom_lengths[chrom]
        chr_bins = int(np.ceil(chr_len / bin_size_bp))

        sub = df[df["chrom"] == chrom]
        if sub.empty:
            continue

        ysum = (sub["brdu_pct"] * sub["Nmod"]).groupby(sub["bin"]).sum()
        wsum = sub["Nmod"].groupby(sub["bin"]).sum()

        mean = ysum / wsum
        mean[wsum < min_cov] = np.nan

        bins_idx = mean.index.to_numpy()
        bins_idx = bins_idx[(bins_idx >= 0) & (bins_idx < chr_bins)]
        H[r, bins_idx] = mean.loc[bins_idx].to_numpy()

    # Covert the precentage into a probability, and quantsize
    H_prob = H / 100.0
    H_prob = np.round(H_prob / 0.25) * 0.25
    H_prob = np.clip(H_prob, 0.0, 1.0)

    return H_prob, x_kb, chroms_used


def plot_heatmap(H, x_kb, chroms_used, out_path: str, title: str):
    fig, ax = plt.subplots(figsize=(16, 6))

    im = ax.imshow(
        H,
        aspect="auto",
        interpolation="nearest",
        origin="upper",
        vmin=0.0,
        vmax=1.0
    )

    ax.set_yticks(np.arange(len(chroms_used)))
    ax.set_yticklabels(chroms_used)

    ax.set_xlabel("Genomic position (kb)")
    ax.set_ylabel("Chromosome")
    ax.set_title(title)

    n = len(x_kb)
    ticks = np.linspace(0, n - 1, 10).astype(int)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{int(x_kb[t])}" for t in ticks])

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("BrdU probability")
    cbar.set_ticks([0.0, 0.25, 0.50, 0.75, 1.0])

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    df = prepare_dataframe()
    output_dir = get_output_dir()

    chroms = sorted(df["chrom"].unique(), key=chrom_sort_key)

    bin_size_bp = 1000
    min_cov = 10.0

    H, x_kb, chroms_used = make_chr_heatmap(
        df,
        chroms=chroms,
        bin_size_bp=bin_size_bp,
        min_cov=min_cov,
        include_extra=False
    )

    out_path = os.path.join(output_dir, f"brdu_heatmap_chr1-16_{bin_size_bp}bp_bins.png")
    plot_heatmap(
        H, x_kb, chroms_used,
        out_path=out_path,
        title=f"M-phase BrdU probability heatmap (chr1–16)"
    )

    print(f"Saved heatmap to: {out_path}")


if __name__ == "__main__":
    main()
