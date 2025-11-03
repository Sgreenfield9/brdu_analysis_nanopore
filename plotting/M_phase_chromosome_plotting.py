import pandas as pd
import matplotlib.pyplot as plt
import os

# -----------------------------
# Load positive and negative bedgraphs
# -----------------------------
pos_df = pd.read_csv(
    "/Users/samuelgreenfiled/Desktop/SG_EXP_2025_010/sort_index/M_phase_pileup/b_positive.bedgraph",
    sep="\t", header=None,
    names=["chrom", "start", "end", "frac_mod", "Nmod"]
)

neg_df = pd.read_csv(
    "/Users/samuelgreenfiled/Desktop/SG_EXP_2025_010/sort_index/M_phase_pileup/b_negative.bedgraph",
    sep="\t", header=None,
    names=["chrom", "start", "end", "frac_mod", "Nmod"]
)

# -----------------------------
# Combine positive and negative
# -----------------------------
df = pd.concat([pos_df, neg_df], ignore_index=True)

# Compute BrdU counts
df["BrdU_count"] = df["frac_mod"] * df["Nmod"]

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
df["chrom"] = df["chrom"].map(genbank_to_chr)
df = df.dropna(subset=["chrom"])
df["chrom"] = df["chrom"].astype(str)

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

# -----------------------------
# Create output folder
# -----------------------------
output_dir = "/Users/samuelgreenfiled/Desktop/SG_EXP_2025_010/sort_index/M_phase_pileup/genome_browser"
os.makedirs(output_dir, exist_ok=True)

# -----------------------------
# Plot each chromosome separately
# -----------------------------
for chrom in sorted(df["chrom"].unique(), key=chrom_sort_key):
    chrom_df = df[df["chrom"] == chrom].copy()
    
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
    ax.legend()
    plt.tight_layout()
    
    # Save figure
    save_path = os.path.join(output_dir, f"chromosome_{chrom}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()

print(f"Plots saved to {output_dir}")