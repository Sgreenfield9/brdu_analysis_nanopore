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
import cairosvg
import numpy as np
from pathlib import Path
from dotenv import load_dotenv
from typing import Optional
from matplotlib.ticker import ScalarFormatter

from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors as rl_colors

# Load positive and negative bedgraphs
REPO_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")


def load_bedgraph_data():
    """
    Load BrdU incorporation data from positive and negative strand bedgraph files.
    """
    pos_path = os.getenv('POSITIVE_BEDGRAPH_M')
    neg_path = os.getenv('NEGATIVE_BEDGRAPH_M')

    if not pos_path or not neg_path:
        raise ValueError("POSITIVE_BEDGRAPH_M and NEGATIVE_BEDGRAPH_M must be set in env/.env")

    pos_df = pd.read_csv(
        pos_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "frac_mod", "Nmod"]
    )

    neg_df = pd.read_csv(
        neg_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "frac_mod", "Nmod"]
    )

    return pos_df, neg_df


def get_output_dir():
    """
    Determine the output directory for saving plots.
    """
    output_dir = os.getenv("OUTPUT_DIR_M")
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "output", "M_phase_pileup", "genome_browser")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def prepare_dataframe():
    """
    Combine the positive and negative strand BrdU data into a single DataFrame. 
    """
    pos_df, neg_df = load_bedgraph_data()
    df = pd.concat([pos_df, neg_df], ignore_index=True)

    df["BrdU_count"] = df["frac_mod"] * df["Nmod"]
    df["chrom"] = df["chrom"].map(map_chromosome)
    df = df.dropna(subset=["chrom"])
    df["chrom"] = df["chrom"].astype(str)
    return df


def load_motifs_from_bed(env_key: str, description: str):
    bed_path = os.getenv(env_key)
    if not bed_path:
        raise ValueError(f"{env_key} must be set in env/.env")
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"{description} BED file not found at {bed_path}")

    motifs_df = pd.read_csv(bed_path, sep="\t", header=None)
    if motifs_df.shape[1] < 3:
        raise ValueError(f"{description} BED file at {bed_path} has fewer than 3 columns")

    motifs_df = motifs_df.iloc[:, :6]
    motifs_df.columns = ["chrom", "start", "end", "name", "score", "strand"][:motifs_df.shape[1]]

    motifs_df["chrom"] = motifs_df["chrom"].map(map_chromosome)
    motifs_df = motifs_df.dropna(subset=["chrom"])
    motifs_df["chrom"] = motifs_df["chrom"].astype(str)
    return motifs_df


def load_g4_motifs():
    g4_path = os.getenv("G4_MOTIFS_BED")
    if g4_path:
        return load_motifs_from_bed("G4_MOTIFS_BED", "G4 motifs")

    g4_fallback = REPO_ROOT / "results" / "G4_extraction" / "g4.motifs.bed"
    if not os.path.exists(g4_fallback):
        raise FileNotFoundError(f"G4 motifs BED file not found at {g4_fallback}")

    g4_df = pd.read_csv(
        g4_fallback,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"]
    )
    g4_df["chrom"] = g4_df["chrom"].map(map_chromosome)
    g4_df = g4_df.dropna(subset=["chrom"])
    g4_df["chrom"] = g4_df["chrom"].astype(str)
    return g4_df


def load_trna_motifs():
    return load_motifs_from_bed("tRNA_MOTIFS_BED", "tRNA motifs")


def load_te_motifs_all():
    return load_motifs_from_bed("TE_MOTIFS_BED_ALL", "TE motifs (all)")


def load_te_motifs_bodies():
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
    if chrom_id in genbank_to_chr:
        return genbank_to_chr[chrom_id]
    if chrom_id in ncbi_refseq_to_chr:
        return ncbi_refseq_to_chr[chrom_id]
    return None


# -----------------------------
# Define subtelomeric regions
# -----------------------------
SUBTEL_SIZE = 30000  # 30 kb from each chromosome end

chrom_lengths = {
    "1": 230218, "2": 813184, "3": 316620, "4": 1531933,
    "5": 576874, "6": 270161, "7": 1090940, "8": 562643,
    "9": 439888, "10": 745751, "11": 666816, "12": 1078177,
    "13": 924431, "14": 784333, "15": 1091291, "16": 948066,
    "MT": 85779, "p2-micron": 6318
}


def chrom_sort_key(x):
    if x.isdigit():
        return (0, int(x))
    else:
        return (1, x)


def smooth_counts(series, window=200):
    return series.rolling(window=window, min_periods=1, center=True).mean()


def trim_horizontal_whitespace(img, thresh=0.985, pad_px=2):
    """
    Trims left/right near-white padding from an RGB/RGBA image.
    Keeps full height. pad_px keeps a small margin so you don't cut ticks.
    """
    if img.ndim != 3:
        return img

    rgb = img[..., :3]
    nonwhite = np.any(rgb < thresh, axis=2)        # not white-ish
    cols = np.where(nonwhite.any(axis=0))[0]       # columns with any drawing

    if cols.size == 0:
        return img

    c0 = max(int(cols[0]) - pad_px, 0)
    c1 = min(int(cols[-1]) + 1 + pad_px, img.shape[1])
    return img[:, c0:c1, :]


def render_feature_diagram_png(
    chrom: str,
    chr_length: int,
    g4_chrom: pd.DataFrame,
    trna_chrom: pd.DataFrame,
    te_chrom: pd.DataFrame,
    out_png: str,
    subtel_size: int = SUBTEL_SIZE,
    width: int = 3000,
    height: int = 260,
):
    """
    Creates a linear GenomeDiagram showing feature tracks for G4/tRNA/TE.

    GenomeDiagram -> SVG
    SVG -> PNG via cairosvg
    Final output is a PNG written to out_png.
    """
    diagram = GenomeDiagram.Diagram(f"Chr {chrom} features")

    # Track 1: ruler + subtelomeres
    track_ruler = diagram.new_track(
        1, name="Ruler", greytrack=False, scale=True, scale_ticks=True
    )
    ruler_set = track_ruler.new_set()

    left = SeqFeature(FeatureLocation(0, min(subtel_size, chr_length)), type="subtelomere")
    ruler_set.add_feature(left, sigil="BOX", color=rl_colors.red, alpha=0.25)

    right_start = max(0, chr_length - subtel_size)
    right = SeqFeature(FeatureLocation(right_start, chr_length), type="subtelomere")
    ruler_set.add_feature(right, sigil="BOX", color=rl_colors.red, alpha=0.25)

    # Track 2: TE
    track_te = diagram.new_track(2, name="TE", greytrack=True)
    te_set = track_te.new_set()
    for _, row in te_chrom.iterrows():
        start, end = int(row["start"]), int(row["end"])
        if end > start:
            feat = SeqFeature(FeatureLocation(start, end), type="TE")
            te_set.add_feature(feat, sigil="BOX", color=rl_colors.orange, alpha=0.8)

    # Track 3: tRNA
    track_trna = diagram.new_track(3, name="tRNA", greytrack=True)
    trna_set = track_trna.new_set()
    for _, row in trna_chrom.iterrows():
        start, end = int(row["start"]), int(row["end"])
        if end > start:
            feat = SeqFeature(FeatureLocation(start, end), type="tRNA")
            trna_set.add_feature(feat, sigil="BOX", color=rl_colors.purple, alpha=0.9)

    # Track 4: G4
    track_g4 = diagram.new_track(4, name="G4", greytrack=True)
    g4_set = track_g4.new_set()
    for _, row in g4_chrom.iterrows():
        start, end = int(row["start"]), int(row["end"])
        if end > start:
            feat = SeqFeature(FeatureLocation(start, end), type="G4")
            g4_set.add_feature(feat, sigil="BOX", color=rl_colors.green, alpha=0.9)

    diagram.draw(
        format="linear",
        pagesize=(width, height),
        start=0,
        end=chr_length,
        fragments=1,
    )

    out_svg = os.path.splitext(out_png)[0] + ".svg"
    diagram.write(out_svg, "SVG")
    cairosvg.svg2png(url=out_svg, write_to=out_png)


def main():
    df = prepare_dataframe()
    g4_df = load_g4_motifs()
    trna_df = load_trna_motifs()
    te_df = load_te_motifs_all()
    # te_df = load_te_motifs_bodies()

    output_dir = get_output_dir()
    tmp_dir = os.path.join(output_dir, "_tmp_feature_diagrams")
    os.makedirs(tmp_dir, exist_ok=True)

    for chrom in sorted(df["chrom"].unique(), key=chrom_sort_key):
        chrom_df = df[df["chrom"] == chrom].copy()
        g4_chrom = g4_df[g4_df["chrom"] == chrom]
        trna_chrom = trna_df[trna_df["chrom"] == chrom]
        te_chrom = te_df[te_df["chrom"] == chrom]

        chrom_df = chrom_df.sort_values("start").drop_duplicates(subset="start")

        chrom_df["BrdU_smooth"] = smooth_counts(chrom_df["BrdU_count"], window=1000)
        chrom_df["Nmod_smooth"] = smooth_counts(chrom_df["Nmod"], window=1000)

        chr_length = int(chrom_lengths.get(chrom, chrom_df["end"].max()))

        # Render GenomeDiagram PNG for this chromosome
        feature_png = os.path.join(tmp_dir, f"chr_{chrom}_features.png")
        render_feature_diagram_png(
            chrom=chrom,
            chr_length=chr_length,
            g4_chrom=g4_chrom,
            trna_chrom=trna_chrom,
            te_chrom=te_chrom,
            out_png=feature_png,
            subtel_size=SUBTEL_SIZE,
        )

        feature_img = plt.imread(feature_png)
        feature_img = trim_horizontal_whitespace(feature_img, thresh=0.985, pad_px=2)

        # Combined figure
        fig = plt.figure(figsize=(15, 6))
        gs = fig.add_gridspec(
            2, 2,
            height_ratios=[3.5, 1.2],
            width_ratios=[0.08, 0.92],
            hspace=0.05,
            wspace=0.0
        )

        # Top plot spans both columns
        ax = fig.add_subplot(gs[0, :])

        ax.axvspan(0, SUBTEL_SIZE, alpha=0.2, color='red', label='Subtelomeric region')
        ax.axvspan(chr_length - SUBTEL_SIZE, chr_length, alpha=0.2, color='red')

        ax.fill_between(
            chrom_df["start"], 0, chrom_df["Nmod_smooth"],
            color='lightgrey', alpha=0.6, label="Coverage (Nmod, smoothed)"
        )
        ax.plot(
            chrom_df["start"], chrom_df["BrdU_smooth"],
            color='blue', linewidth=1, label="BrdU count (smoothed)"
        )

        ax.set_ylabel("Read count")
        ax.set_title(f"Mitosis BrdU pileup along chromosome {chrom} (subtelomeres highlighted)")
        ax.set_xlim(0, chr_length)
        ax.set_ylim(0, 22)
        ax.set_xlabel("Genomic position (bp)")

        # Make x-axis read as basepairs (no 1e6 offset)
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.ticklabel_format(style='plain', axis='x')

        ax.grid(True, axis="x", alpha=0.15)
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)

        # Bottom-left: labels
        ax_lab = fig.add_subplot(gs[1, 0])
        ax_lab.axis("off")
        ax_lab.text(0.98, 0.80, "G4",    ha="right", va="center", fontsize=11, transform=ax_lab.transAxes)
        ax_lab.text(0.98, 0.55, "tRNA",  ha="right", va="center", fontsize=11, transform=ax_lab.transAxes)
        ax_lab.text(0.98, 0.30, "TE",    ha="right", va="center", fontsize=11, transform=ax_lab.transAxes)
        ax_lab.text(0.98, 0.08, "Subtel", ha="right", va="center", fontsize=10, transform=ax_lab.transAxes)

        # Bottom-right: GenomeDiagram image aligned to genomic coordinates
        ax_img = fig.add_subplot(gs[1, 1], sharex=ax)
        ax_img.imshow(
            feature_img,
            aspect="auto",
            extent=[0, chr_length, 0, 1],
            interpolation="nearest"
        )
        ax_img.set_xlim(0, chr_length)
        ax_img.axis("off")
        ax_img.text(
            0.5, 1.02, "Genome features (GenomeDiagram)",
            ha="center", va="bottom",
            transform=ax_img.transAxes,
            fontsize=10
        )

        plt.tight_layout()

        save_path = os.path.join(output_dir, f"chromosome_{chrom}.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()

    print(f"Plots saved to {output_dir}")


if __name__ == "__main__":
    main()
