"""
This script visualizes BrdU incorporation acroos yeast chromosomes
during M phase (mitosis). 

The script:
1. Loads BrdU bedgraph data
2. Overlays genomic features (G4, motifs, tRNA's, transposable elements (TE))
3. Highlights subtelomeric regions
4. Generate per-chromosome browser plots 
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
from dotenv import load_dotenv
from typing import Optional

# Load positive and negative bedgraphs
REPO_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")

def load_bedgraph_data():
    """
    Load BrdU incorporation data from positive and negative strand bedgraph files.
    """
    # Get paths
    pos_path = os.getenv('POSITIVE_BEDGRAPH_M')
    neg_path = os.getenv('NEGATIVE_BEDGRAPH_M')

    if not pos_path or not neg_path:
        raise ValueError("POSITIVE_BEDGRAPH_M and NEGATIVE_BEDGRAPH_M must be set in env/.env")
    
    # Read postivie strand data
    pos_df = pd.read_csv(
        pos_path,
        sep="\t",
        header = None,
        names = ["chrom", "start", "end", "frac_mod", "Nmod"]
    )

    # Read negative strand data
    neg_df = pd.read_csv(
        neg_path,
        sep="\t",
        header = None,
        names = ["chrom", "start", "end", "frac_mod", "Nmod"]
    )

    return pos_df, neg_df

def get_output_dir():
    """
    Determine the output directory for saving plots.

    Checks for OUTPUT_DIR_M environment variable. If not set, creates a default
    directory structure: output/M_phase_pileup/genome_browser/ 
    """
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
    """
    Combine the positive and negative strand BrdU data into a single DataFrame. 
    """
    pos_df, neg_df = load_bedgraph_data()
    df = pd.concat([pos_df, neg_df], ignore_index=True)
    
    # Compute BrdU counts
    df["BrdU_count"] = df["frac_mod"] * df["Nmod"]
    df["chrom"] = df["chrom"].map(map_chromosome)
    df = df.dropna(subset=["chrom"])
    df["chrom"] = df["chrom"].astype(str)
    return df

def load_motifs_from_bed(env_key: str, description: str):
    """
    Load genomic feature coordinates from a BED file; this is specified in the environment variables.
    """
    # Get BED file path from environment variable 
    bed_path = os.getenv(env_key)
    if not bed_path:
        raise ValueError(f"{env_key} must be set in env/.env")
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"{description} BED file not found at {bed_path}")

    # Read BED file
    motifs_df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None
    )
    # Make sure BED has three columns (chrom, start, end)
    if motifs_df.shape[1] < 3:
        raise ValueError(f"{description} BED file at {bed_path} has fewer than 3 columns")
    # Some BEDs include extra columns (e.g., feature type). Keep the first 6 for plotting.
    motifs_df = motifs_df.iloc[:, :6]
    motifs_df.columns = ["chrom", "start", "end", "name", "score", "strand"][:motifs_df.shape[1]]
    motifs_df["chrom"] = motifs_df["chrom"].map(map_chromosome)
    motifs_df = motifs_df.dropna(subset=["chrom"])
    motifs_df["chrom"] = motifs_df["chrom"].astype(str)
    return motifs_df

def load_g4_motifs():
    """
    Load G4 quadruplex motif coordinates from BED file. Returns
    a DataFrame contatin G4 motif coordinates with standardized columns
    """
    # Load from the environment variable path
    g4_path = os.getenv("G4_MOTIFS_BED")
    if g4_path:
        return load_motifs_from_bed("G4_MOTIFS_BED", "G4 motifs")

    g4_fallback = REPO_ROOT / "results" / "G4_extraction" / "g4.motifs.bed"
    if not os.path.exists(g4_fallback):
        raise FileNotFoundError(f"G4 motifs BED file not found at {g4_fallback}")

    # Read fallback fule with explicit column names
    g4_df = pd.read_csv(
        g4_fallback,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"]
    )

    # Map and filter chromosomes
    g4_df["chrom"] = g4_df["chrom"].map(map_chromosome)
    g4_df = g4_df.dropna(subset=["chrom"])
    g4_df["chrom"] = g4_df["chrom"].astype(str)
    return g4_df

def load_trna_motifs():
    """
    Load tRNA (transfer RNA) feature coordinates from BED file. 
    """
    return load_motifs_from_bed("tRNA_MOTIFS_BED", "tRNA motifs")

def load_te_motifs_all():
    """
    Load all transposable elements (TE) features including bodies and LTRs.
    """
    return load_motifs_from_bed("TE_MOTIFS_BED_ALL", "TE motifs (all)")

def load_te_motifs_bodies():
    """
    Load transposable element bodies only (excluding LTRs)
    """
    return load_motifs_from_bed("TE_MOTIFS_BED_BODIES", "TE motifs (bodies)")

# -----------------------------
# Map GenBank IDs to chromosome numbers
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

# -----------------------------
# Define subtelomeric regions (in bp from each end)
# You can adjust this value based on your organism
# For S. cerevisiae, subtelomeres are typically 0-30kb from each end
# -----------------------------
SUBTEL_SIZE = 30000  # 30 kb from each chromosome end

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

def main():
    df = prepare_dataframe()
    g4_df = load_g4_motifs()
    trna_df = load_trna_motifs()
    te_df = load_te_motifs_all()
    # te_df = load_te_motifs_bodies()

    output_dir = get_output_dir()

    # -----------------------------
    # Plot each chromosome separately
    # -----------------------------
    for chrom in sorted(df["chrom"].unique(), key=chrom_sort_key):
        chrom_df = df[df["chrom"] == chrom].copy()
        g4_chrom = g4_df[g4_df["chrom"] == chrom]
        trna_chrom = trna_df[trna_df["chrom"] == chrom]
        te_chrom = te_df[te_df["chrom"] == chrom]
        
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
        ax.set_xlim(0, chr_length)

        # Extend y-limit to make room for motif tracks above the data
        ax.set_ylim(0, 22)
        y_max = 20  # Keep motifs relative to data range, not extended limit

        # Position motif tracks above the main data
        g4_line_height = y_max * 1.10
        trna_line_height = y_max * 1.05  
        te_line_height = y_max * 1.00
        band_height = y_max * 0.04  # Thicker bands for visibility

        # Plot G4 motifs as thick horizontal bands
        first_g4 = True
        for _, motif in g4_chrom.iterrows():
            label = "G4 motif" if first_g4 else None
            ax.axvspan(
                motif["start"],
                motif["end"],
                ymin=(g4_line_height - band_height/2) / 22,
                ymax=(g4_line_height + band_height/2) / 22,
                color="green",
                alpha=0.8,
                label=label
            )
            first_g4 = False

        # Plot tRNA motifs as thick horizontal bands
        first_trna = True
        for _, motif in trna_chrom.iterrows():
            label = "tRNA motif" if first_trna else None
            ax.axvspan(
                motif["start"],
                motif["end"],
                ymin=(trna_line_height - band_height/2) / 22,
                ymax=(trna_line_height + band_height/2) / 22,
                color="purple",
                alpha=0.8,
                label=label
            )
            first_trna = False

        # Plot TE motifs as thick horizontal bands
        first_te = True
        for _, motif in te_chrom.iterrows():
            label = "TE motif" if first_te else None
            ax.axvspan(
                motif["start"],
                motif["end"],
                ymin=(te_line_height - band_height/2) / 22,
                ymax=(te_line_height + band_height/2) / 22,
                color="orange",
                alpha=0.8,
                label=label
            )
            first_te = False

        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
        plt.tight_layout()
        
        # Save figure
        save_path = os.path.join(output_dir, f"chromosome_{chrom}.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()

    print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    main()
