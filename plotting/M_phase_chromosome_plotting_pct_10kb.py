import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from dotenv import load_dotenv
from matplotlib.ticker import ScalarFormatter

import M_phase_chromosome_plotting as m

REPO_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")


def get_output_dir():
    """
    This gets the output directory for where we want to put our new plots. This will
    be given in the env/e.env
    """
    output_dir = os.getenv("BRDU_M_PCT_10KB_OUTPUT")
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "output", "M_Phase_pct_10kb")
    os.makedirs(output_dir, exist_ok=True)
    return Path(output_dir)


def load_M_pct_10kb():
    """
    This method finds the path for the dataset we are using. If the dataset is not found
    it warns the user to put it in the env/.env. If the path does not exist it will warn
    the user that the path does not exist.
    """
    M_10kb_pct = os.getenv("BRDU_M_PCT_10KB")
    if not M_10kb_pct:
        raise ValueError("BRDU_M_PCT_10KB must be set in the env/.env")
    if not os.path.exists(M_10kb_pct):
        raise FileNotFoundError(f"BRDU_M_PCT_10KB not found at {M_10kb_pct}")
    return pd.read_csv(M_10kb_pct)


# Input come fron env/.env
poi = load_M_pct_10kb().head(20)
out_dir = get_output_dir()

# Getting our data from the M_phase_chromosome_plotting.py
df = m.prepare_dataframe()
g4_df = m.load_g4_motifs()
trna_df = m.load_trna_motifs()
te_df = m.load_te_motifs_all()

# Cache per-chromosome smoothed dataframe once
chrom_cache = {}  # chrom -> (chrom_df_smoothed, chr_length)

for _, row in poi.iterrows():
    chrom = str(row["chrom"])
    center = int(row["start"])

    # Build / reuse cached chromosome dataframe (sorted + smoothed) once per chrom
    if chrom not in chrom_cache:
        chrom_df = df[df["chrom"] == chrom].copy()
        if chrom_df.empty:
            chrom_cache[chrom] = (None, None)
        else:
            chrom_df = chrom_df.sort_values("start").drop_duplicates(subset="start")

            # Smooth + compute % ONCE per chromosome
            chrom_df["BrdU_smooth"] = m.smooth_counts(chrom_df["BrdU_count"], window=1000)
            chrom_df["Nmod_smooth"] = m.smooth_counts(chrom_df["Nmod"], window=1000)
            chrom_df["BrdU_pct"] = np.where(
                chrom_df["Nmod"] > 0,
                100.0 * chrom_df["BrdU_count"] / chrom_df["Nmod"],
                np.nan
            )
            chrom_df.loc[chrom_df["BrdU_pct"] <= 0, "BrdU_pct"] = np.nan
            chrom_df["BrdU_pct_smooth"] = m.smooth_counts(chrom_df["BrdU_pct"], window=1000)

            chr_length = int(m.chrom_lengths.get(chrom, chrom_df["end"].max()))
            chrom_cache[chrom] = (chrom_df, chr_length)

    chrom_df, chr_length = chrom_cache.get(chrom, (None, None))
    if chrom_df is None or chr_length is None:
        continue

    # Window boundaries: 10kb total (5kb each side), with left-edge handling
    if center < 5000:
        left = 0
        right = min(chr_length, center + 10000)
    else:
        left = max(0, center - 5000)
        right = min(chr_length, center + 5000)

    # Slice to just the window (fast)
    window_mask = (chrom_df["start"] >= left) & (chrom_df["start"] <= right)
    chrom_window = chrom_df.loc[window_mask].copy()
    if chrom_window.empty:
        continue

    # Slice feature tracks to the window
    g4_chrom = g4_df[(g4_df["chrom"] == chrom) & (g4_df["end"] >= left) & (g4_df["start"] <= right)]
    trna_chrom = trna_df[(trna_df["chrom"] == chrom) & (trna_df["end"] >= left) & (trna_df["start"] <= right)]
    te_chrom = te_df[(te_df["chrom"] == chrom) & (te_df["end"] >= left) & (te_df["start"] <= right)]

    # Layout
    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(
        3, 1,
        height_ratios=[3.5, 0.40, 2.5],
        hspace=0.15
    )

    ax = fig.add_subplot(gs[0, 0])

    # Subtel background (still drawn using full-chr coordinates; xlim will clip)
    ax.axvspan(0, m.SUBTEL_SIZE, alpha=0.2, color='red', label='Subtelomeric region')
    ax.axvspan(chr_length - m.SUBTEL_SIZE, chr_length, alpha=0.2, color='red')

    # Plot ONLY the window
    ax.fill_between(
        chrom_window["start"], 0, chrom_window["Nmod_smooth"],
        color='lightgrey', alpha=0.6, label="Coverage (Nmod, smoothed)"
    )
    ax.plot(
        chrom_window["start"], chrom_window["BrdU_smooth"],
        color='blue', linewidth=1, label="BrdU count (smoothed)"
    )

    ax.set_ylabel("Read count")
    ax.set_title(f"Mitosis BrdU pileup chr{chrom} window {left}-{right} (center {center})")
    ax.set_ylim(0, 22)
    ax.set_xlabel("Genomic position (bp)")
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax.ticklabel_format(style='plain', axis='x')
    ax.grid(True, axis="x", alpha=0.15)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)

    # SPEED + LOOK: compute bin_size from WINDOW, and keep it sane
    window_len = max(1, right - left)

    # This controls how "continuous" the stick rulers look
    approx_sticks = 2000

    # Never let bins get too tiny (1bp bins = many vlines = slow)
    bin_size = max(10, int(window_len // approx_sticks))  # 10bp minimum
    bin_size = min(bin_size, 100)                         # 100bp maximum (for safety)

    m.add_all_rulers_panel(
        fig=fig,
        gs_slot=gs[2, 0],
        ax_sharex=ax,
        chrom_df=chrom_window,  # use window df for rulers
        chr_length=chr_length,
        g4_chrom=g4_chrom,
        trna_chrom=trna_chrom,
        te_chrom=te_chrom,
        bin_size=bin_size
    )

    # Force window x-limits on all axes
    for a in fig.axes:
        a.set_xlim(left, right)

    save_path = out_dir / f"chromosome_{chrom}_start_{center}_10kb.png"

    # Faster saves while iterating
    fig.subplots_adjust(right=0.80)
    fig.savefig(save_path, dpi=250)  
    plt.close(fig)

print("wrote", out_dir)
