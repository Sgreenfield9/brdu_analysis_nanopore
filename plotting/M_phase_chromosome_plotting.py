import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
from dotenv import load_dotenv
from typing import Optional

# -----------------------------
# Load positive and negative bedgraphs
# -----------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")

# -----------------------------
# Chromosome mapping helpers
# -----------------------------
genbank_to_chr = {
    "CM007964.1": "1",  "CM007965.1": "2",  "CM007966.1": "3",  "CM007967.1": "4",
    "CM007968.1": "5",  "CM007969.1": "6",  "CM007970.1": "7",  "CM007971.1": "8",
    "CM007972.1": "9",  "CM007973.1": "10", "CM007974.1": "11", "CM007975.1": "12",
    "CM007976.1": "13", "CM007977.1": "14", "CM007978.1": "15", "CM007979.1": "16",
    "CM007980.1": "p2-micron", "CM007981.1": "MT"
}

ncbi_refseq_to_chr = {
    "NC_001133.9": "1",   "NC_001134.8": "2",   "NC_001135.5": "3",   "NC_001136.10": "4",
    "NC_001137.3": "5",   "NC_001138.5": "6",   "NC_001139.9": "7",   "NC_001140.6": "8",
    "NC_001141.2": "9",   "NC_001142.9": "10",  "NC_001143.9": "11",  "NC_001144.5": "12",
    "NC_001145.3": "13",  "NC_001146.8": "14",  "NC_001147.6": "15",  "NC_001148.4": "16",
    "NC_001224.1": "MT"
}

def map_chromosome(chrom_id: str) -> Optional[str]:
    """Map assembly-specific chromosome IDs to numeric/short labels used in plots."""
    if chrom_id in genbank_to_chr:
        return genbank_to_chr[chrom_id]
    if chrom_id in ncbi_refseq_to_chr:
        return ncbi_refseq_to_chr[chrom_id]
    return None

def load_bedgraph_data():

    # Get paths
    pos_path = os.getenv('POSITIVE_BEDGRAPH_M')
    neg_path = os.getenv('NEGATIVE_BEDGRAPH_M')

    if not pos_path or not neg_path:
        raise ValueError("POSITIVE_BEDGRAPH_M and NEGATIVE_BEDGRAPH_M must be set in env/.env")
    
    # Read Data
    pos_df = pd.read_csv(
        pos_path,
        sep="\t",
        header = None,
        names = ["chrom", "start", "end", "frac_mod", "Nmod"]
    )

    neg_df = pd.read_csv(
        neg_path,
        sep="\t",
        header = None,
        names = ["chrom", "start", "end", "frac_mod", "Nmod"]
    )

    return pos_df, neg_df

def get_output_dir():
    # Allow overriding via env; otherwise use repo-local output folder
    output_dir = os.getenv("OUTPUT_DIR_M")
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "output", "M_phase_pileup", "genome_browser")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

# -----------------------------
# Combine positive and negative
# -----------------------------
def prepare_dataframe():
    pos_df, neg_df = load_bedgraph_data()
    df = pd.concat([pos_df, neg_df], ignore_index=True)
    
    # Compute BrdU counts
    df["BrdU_count"] = df["frac_mod"] * df["Nmod"]
    df["chrom"] = df["chrom"].map(map_chromosome)
    df = df.dropna(subset=["chrom"])
    df["chrom"] = df["chrom"].astype(str)
    return df

# -----------------------------
# Load G4 motifs
# -----------------------------
def load_g4_motifs():
    g4_path = os.getenv("G4_MOTIFS_BED")
    if not g4_path:
        g4_path = REPO_ROOT / "results" / "G4_extraction" / "g4.motifs.bed"
    if not os.path.exists(g4_path):
        raise FileNotFoundError(f"G4 motifs BED file not found at {g4_path}")

    g4_df = pd.read_csv(
        g4_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"]
    )
    g4_df["chrom"] = g4_df["chrom"].map(map_chromosome)
    g4_df = g4_df.dropna(subset=["chrom"])
    g4_df["chrom"] = g4_df["chrom"].astype(str)
    return g4_df

# -----------------------------
# Define subtelomeric regions (in bp from each end)
# You can adjust this value based on your organism
# For S. cerevisiae, subtelomeres are typically 0-30kb from each end
# -----------------------------
SUBTEL_SIZE = 30000  # 30 kb from each chromosome end
DEFAULT_G4_ZOOM_PADDING_BP = 2000  # 2 kb window around each G4 motif for zoom plots

# Chromosome lengths for S. cerevisiae (adjust if using different organism)
chrom_lengths = {
    "1": 230218, "2": 813184, "3": 316620, "4": 1531933,
    "5": 576874, "6": 270161, "7": 1090940, "8": 562643,
    "9": 439888, "10": 745751, "11": 666816, "12": 1078177,
    "13": 924431, "14": 784333, "15": 1091291, "16": 948066,
    "MT": 85779, "p2-micron": 6318
}

# Sorting key for chromosomes
def chrom_sort_key(x):
    if x.isdigit():
        return (0, int(x))
    else:
        return (1, x)

# -----------------------------
# Smoothing function
# -----------------------------
def smooth_counts(series, window=200):
    return series.rolling(window=window, min_periods=1, center=True).mean()

def get_g4_zoom_padding_bp():
    value = os.getenv("G4_ZOOM_PADDING_BP", str(DEFAULT_G4_ZOOM_PADDING_BP))
    try:
        return int(value)
    except ValueError:
        return DEFAULT_G4_ZOOM_PADDING_BP

def main():
    df = prepare_dataframe()
    g4_df = load_g4_motifs()

    output_dir = get_output_dir()

    # -----------------------------
    # Plot each chromosome separately
    # -----------------------------
    for chrom in sorted(df["chrom"].unique(), key=chrom_sort_key):
        chrom_df = df[df["chrom"] == chrom].copy()
        g4_chrom = g4_df[g4_df["chrom"] == chrom]
        
        # Sort by position
        chrom_df = chrom_df.sort_values("start")
        
        # Optional: remove duplicates
        chrom_df = chrom_df.drop_duplicates(subset="start")
        
        # Smooth per chromosome
        chrom_df["BrdU_smooth"] = smooth_counts(chrom_df["BrdU_count"], window=1000)
        chrom_df["Nmod_smooth"] = smooth_counts(chrom_df["Nmod"], window=1000)
        
        # Get chromosome length
        chr_length = chrom_lengths.get(chrom, chrom_df["end"].max())
        
        # Create plot
        fig, ax = plt.subplots(figsize=(15, 4))
        
        # Highlight subtelomeric regions
        # Left subtelomere (0 to SUBTEL_SIZE)
        ax.axvspan(0, SUBTEL_SIZE, alpha=0.2, color='red', label='Subtelomeric region')
        
        # Right subtelomere (chr_length - SUBTEL_SIZE to chr_length)
        ax.axvspan(chr_length - SUBTEL_SIZE, chr_length, alpha=0.2, color='red')
        
        # Plot data
        ax.fill_between(chrom_df["start"], 0, chrom_df["Nmod_smooth"],
                        color='lightgrey', alpha=0.6, label="Coverage (Nmod, smoothed)")
        ax.plot(chrom_df["start"], chrom_df["BrdU_smooth"],
                color='blue', linewidth=1, label="BrdU count (smoothed)")

        ax.set_xlabel("Genomic position (bp)")
        ax.set_ylabel("Read count")
        ax.set_title(f"Mitosis BrdU pileup along chromosome {chrom} (subtelomeres highlighted)")
        ax.set_ylim(0, 20)
        ax.set_xlim(0, chr_length)

        # Overlay G4 motifs as green segments across their start/end span.
        # Compute the line height after setting axis limits so segments stay in view.
        g4_line_height = ax.get_ylim()[1] * 0.95
        first_g4 = True
        for _, motif in g4_chrom.iterrows():
            label = "G4 motif" if first_g4 else None
            ax.plot(
                [motif["start"], motif["end"]],
                [g4_line_height, g4_line_height],
                color="green",
                linewidth=2.5,
                solid_capstyle="butt",
                label=label
            )
            first_g4 = False
        ax.legend()
        plt.tight_layout()
        
        # Save figure
        save_path = os.path.join(output_dir, f"chromosome_{chrom}.png")
        plt.savefig(save_path, dpi=300)
        plt.close()

        # -----------------------------
        # Zoomed plots per G4 motif
        # -----------------------------
        if not g4_chrom.empty:
            zoom_padding = get_g4_zoom_padding_bp()
            zoom_dir = os.path.join(output_dir, "g4_zoom")
            os.makedirs(zoom_dir, exist_ok=True)
            for idx, motif in g4_chrom.reset_index(drop=True).iterrows():
                motif_start = int(motif["start"])
                motif_end = int(motif["end"])
                zoom_start = max(0, motif_start - zoom_padding)
                zoom_end = min(chr_length, motif_end + zoom_padding)

                window_df = chrom_df[
                    (chrom_df["start"] >= zoom_start) & (chrom_df["start"] <= zoom_end)
                ]
                if window_df.empty:
                    continue

                fig, ax = plt.subplots(figsize=(10, 4))

                ax.fill_between(
                    window_df["start"],
                    0,
                    window_df["Nmod_smooth"],
                    color="lightgrey",
                    alpha=0.6,
                    label="Coverage (Nmod, smoothed)"
                )
                ax.plot(
                    window_df["start"],
                    window_df["BrdU_smooth"],
                    color="blue",
                    linewidth=1,
                    label="BrdU count (smoothed)"
                )

                local_ymax = max(
                    window_df["Nmod_smooth"].max(),
                    window_df["BrdU_smooth"].max()
                )
                if pd.isna(local_ymax) or local_ymax <= 0:
                    local_ymax = 1
                ax.set_ylim(0, local_ymax * 1.1)
                ax.set_xlim(zoom_start, zoom_end)
                ax.set_xlabel("Genomic position (bp)")
                ax.set_ylabel("Read count")
                ax.set_title(
                    f"G4 zoom for chromosome {chrom}: {motif_start}-{motif_end} bp"
                )

                g4_line_height = ax.get_ylim()[1] * 0.95
                ax.plot(
                    [motif_start, motif_end],
                    [g4_line_height, g4_line_height],
                    color="green",
                    linewidth=3,
                    solid_capstyle="butt",
                    label="G4 motif"
                )
                ax.axvline(motif_start, color="green", linestyle="--", linewidth=1)
                ax.axvline(motif_end, color="green", linestyle="--", linewidth=1)
                ax.text(
                    motif_start,
                    g4_line_height,
                    f"{motif_start}",
                    rotation=90,
                    va="bottom",
                    ha="right",
                    fontsize=8
                )
                ax.text(
                    motif_end,
                    g4_line_height,
                    f"{motif_end}",
                    rotation=90,
                    va="bottom",
                    ha="left",
                    fontsize=8
                )

                ax.legend()
                plt.tight_layout()

                zoom_path = os.path.join(
                    zoom_dir,
                    f"chromosome_{chrom}_g4_{motif_start}_{motif_end}.png"
                )
                plt.savefig(zoom_path, dpi=300)
                plt.close()

    print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    main()