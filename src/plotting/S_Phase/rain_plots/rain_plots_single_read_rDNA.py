import os
from pathlib import Path
from dotenv import load_dotenv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 

REPO_ROOT = Path(__file__).resolve().parents[4]
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")


def load_rain_plot_input():
    """
    Finds our rainplot input file path. It will then convert
    that file into a pandas dataframe to make it easier/faster to plot.
    """
    rain_plot_path = os.getenv("INPUT_FILE_RP_rDNA")

    if not rain_plot_path:
        raise ValueError("INPUT_FILE_RP_rDNA must be set in the env/.env")

    rain_plot_df = pd.read_csv(
        rain_plot_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "read_id", "mod_base", "mod_prob"],
    )

    return rain_plot_df


def load_rfb_coords():
    """
    Load RFB coordinates per read. Expected columns:
    read_id, chrom, strand, read_pos_start, read_pos_end, ref_start, ref_end
    """
    rfb_path = REPO_ROOT / "data" / "S_Phase" / "rDNA" / "fasta" / "RFB_coords.tsv"
    if not rfb_path.exists():
        return pd.DataFrame()

    rfb_df = pd.read_csv(rfb_path, sep="\t")
    if "read_id" not in rfb_df.columns:
        return pd.DataFrame()

    # Normalize and coerce numeric columns
    rfb_df["read_id"] = rfb_df["read_id"].astype(str)
    for col in ["read_pos_start", "read_pos_end", "ref_start", "ref_end"]:
        if col in rfb_df.columns:
            rfb_df[col] = pd.to_numeric(rfb_df[col], errors="coerce")

    return rfb_df


def get_output_dir():
    """
    Determine the output directory for the rain plots
    """
    output_dir = os.getenv("OUTPUT_FILE_RP_rDNA")

    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "output", "rain_plot")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


# Map GenBank IDs to chromosome numbers
genbank_to_chr = {
    "CM007964.1": "1",  "CM007965.1": "2",  "CM007966.1": "3",  "CM007967.1": "4",
    "CM007968.1": "5",  "CM007969.1": "6",  "CM007970.1": "7",  "CM007971.1": "8",
    "CM007972.1": "9",  "CM007973.1": "10", "CM007974.1": "11", "CM007975.1": "12",
    "CM007976.1": "13", "CM007977.1": "14", "CM007978.1": "15", "CM007979.1": "16",
    "CM007980.1": "p2-micron", "CM007981.1": "MT"
}

def plot_rainplots_per_read():
    """
    Rain plots (up to first 100 reads) using adaptive/equal count binning.

    Idea:
      - Instead of fixed-width bins (e.g., every 0.2 kb), we create bins that each contain
        roughly the same number of points (bases) from the read.
      - Dense regions -> narrower bins (more x-resolution)
      - Sparse regions -> wider bins (more smoothing)

    Plot:
      - Optional scatter at bin centers (black)
      - Step curve using ax.stairs(values, edges, fill=False)
    """
    # Using input function to get the df for plotting
    df = load_rain_plot_input()
    rfb_df = load_rfb_coords()
    # Stating where our output directory is located
    outdir = get_output_dir()

    # memory optimization + cleanup
    df["start"] = pd.to_numeric(df["start"], errors="coerce", downcast="integer")
    df["mod_prob"] = pd.to_numeric(df["mod_prob"], errors="coerce", downcast="float")
    df = df.dropna(subset=["read_id", "start", "mod_prob"])

    # Dataset is limited to 100 reads, but in the future we can limit as well
    # We can leave this line of code and comment it out when needed
    # We can also change the value to whatever we'd like
    max_reads = 100

    
    min_T_count = 800 # change this number to whatever threshold you want

    # This is the code that makes it random each time you run the script.
    # Comment this block out if you want the first 100 reads every time.
    rng = np.random.default_rng()

    # Normalize mod_base in case of lowercase / weird formatting
    df["mod_base"] = df["mod_base"].astype(str).str.upper()

    # First: compute T counts per read, then filter to reads that pass threshold
    per_read = (
        df.groupby("read_id", as_index=False)
          .agg(
              T_count=("mod_base", lambda s: (s == "T").sum())
          )
    )

    eligible_ids = per_read.loc[per_read["T_count"] >= min_T_count, "read_id"].to_numpy()

    if eligible_ids.size == 0:
        raise ValueError(
            f"No reads passed the T threshold. min_T_count={min_T_count}. "
            f"Lower min_T_count or verify mod_base contains 'T'."
        )

    # Then: filter those eligible reads so we do not get flat lines
    points_per_bin = 25  # Can change to 100 for more detail, 500 for smoother plots
    min_q_spread = 0.20 # was 0.30
    min_bin_range = 0.12 # as 0.20

    df_eligible = df[df["read_id"].isin(eligible_ids)]
    var_keep_ids = []

    for rid, sub in df_eligible.groupby("read_id", sort=False):
        sub = sub.sort_values("start", kind="mergesort")
        y = sub["mod_prob"].to_numpy(dtype=float)

        if y.size < 10:
            continue

        q05 = float(np.quantile(y, 0.05))
        q95 = float(np.quantile(y, 0.95))
        q_spread = q95 - q05
        if q_spread < min_q_spread:
            continue

        n = y.size
        n_bins = max(1, int(np.ceil(n / points_per_bin)))
        idx_chunks = np.array_split(np.arange(n), n_bins)
        y_vals_tmp = []
        for idx in idx_chunks:
            if idx.size == 0:
                continue
            y_vals_tmp.append(float(np.mean(y[idx])))

        if len(y_vals_tmp) < 2:
            continue

        bin_range = max(y_vals_tmp) - min(y_vals_tmp)
        if bin_range < min_bin_range:
            continue

        var_keep_ids.append(rid)

    eligible_ids = np.asarray(var_keep_ids, dtype=object)

    if eligible_ids.size == 0:
        raise ValueError(
            f"No reads passed the variability filter. "
            f"Try lowering min_q_spread={min_q_spread} or min_bin_range={min_bin_range}."
        )

    # Restrict to reads that have RFB coordinates
    if rfb_df.empty:
        raise ValueError("RFB_coords.tsv is missing or empty; cannot plot RFB lines.")

    rfb_read_ids = set(rfb_df["read_id"].astype(str))
    eligible_ids = np.asarray([rid for rid in eligible_ids if str(rid) in rfb_read_ids], dtype=object)

    if eligible_ids.size == 0:
        raise ValueError(
            "No reads passed filters AND had RFB coords. "
            "Lower thresholds or verify RFB_coords.tsv read_ids match the rain input."
        )

    # Then: sample up to max_reads from those eligible reads
    if eligible_ids.size > max_reads:
        sampled_ids = rng.choice(eligible_ids, size=max_reads, replace=False)
    else:
        sampled_ids = eligible_ids

    df = df[df["read_id"].isin(sampled_ids)]

    # Binning control
    show_scatter = True   # Show the scatter, we can change to False if we don't want scatter

   # Iterate through each read (we have 100 reads in our dataset)
   # Each read is given a "group id", we start at 1 and go to 100
   # This helps us sort the reads, along with making sure each read
   # is unique and has no duplicates
    for i, (_, sub) in enumerate(df.groupby("read_id", sort=False), start=1):
        if i > max_reads:
            break

        # Sort bases within each read by genomic start position
        # to ensure it increases from left-to-right
        # This helps with ordering
        sub = sub.sort_values("start", kind="mergesort")

        # Assinging the genbank ids to a chr_label
        # Making sure there are no duplicates and that all 
        # chromosomes are present
        genbank_ids = sub["chrom"].astype(str).unique()
        if len(genbank_ids) == 1:
            chr_label = genbank_to_chr.get(genbank_ids[0], genbank_ids[0])
        else:
            chr_label = "mixed/" + ",".join(genbank_to_chr.get(g, g) for g in genbank_ids)

        # X axis: kb position within this read (relative)
        read_start = sub["start"].min()
        x = (sub["start"].to_numpy() - read_start) / 1000.0
        y = sub["mod_prob"].to_numpy()

        # Get the number of bases in a genomic read
        # If less than two we will skip 
        n = len(x)
        if n < 2:
            continue

        # Number of bins (at least 1)
        n_bins = max(1, int(np.ceil(n / points_per_bin)))

        # Split indices into equal-count chunks
        idx_chunks = np.array_split(np.arange(n), n_bins)

        # Build adaptive edges and per-bin y values
        edges = [x[0]]
        y_vals = []
        centers = []

        for idx in idx_chunks:
            if idx.size == 0:
                continue

            # ensure chunk x is monotonic (it is, since sub sorted by start)
            # What I mean by monotonic is that it's not really increasing
            # or decreasing (staying consistent)
            x_left = x[idx[0]]
            x_right = x[idx[-1]]

            # If multiple points share the same x, x_right may equal x_left; that's okay,
            # but stairs prefers non-decreasing edges.
            edges.append(x_right)

            # Representative y for the bin (Currently mean, can change to median if needed)
            y_bin = float(np.mean(y[idx]))
            y_vals.append(y_bin)

            # Bin center for optional scatter
            centers.append((x_left + x_right) / 2.0)

        # Convert to numpy
        edges = np.asarray(edges, dtype=float)
        y_vals = np.asarray(y_vals, dtype=float)
        centers = np.asarray(centers, dtype=float)

        # Can produce duplicates, so we enforce this safely
        if len(edges) < 2 or len(y_vals) < 1:
            continue

        # If edges length doesn't match y_vals+1, rebuild edges from chunk boundaries precisely
        if len(edges) != len(y_vals) + 1:
            edges = [x[idx_chunks[0][0]]]
            for idx in idx_chunks:
                if idx.size == 0:
                    continue
                edges.append(x[idx[-1]])
            edges = np.asarray(edges, dtype=float)

        # If we still have mismatch (rare), skip safely
        if len(edges) != len(y_vals) + 1:
            continue

        # Ensure non-decreasing edges
        edges = np.maximum.accumulate(edges)

        # Avoid a completely flat final edge equal to previous
        # Add a tiny epsilon so the last step is drawable
        if edges[-1] == edges[-2]:
            edges[-1] = edges[-1] + 1e-9

        fig, ax = plt.subplots(figsize=(16, 4))

        # Optional scatter
        if show_scatter:
            ax.scatter(x, y, s=1, color="black", alpha=0.3)

        # RFB region (per read) as dashed red vertical lines
        if not rfb_df.empty:
            rid = str(sub["read_id"].iloc[0])
            rfb_row = rfb_df[rfb_df["read_id"] == rid]
            if not rfb_row.empty:
                rfb_row = rfb_row.iloc[0]
                rfb_start = rfb_row.get("read_pos_start", np.nan)
                rfb_end = rfb_row.get("read_pos_end", np.nan)

                # Fallback to reference coords if read positions are missing
                if pd.isna(rfb_start) or pd.isna(rfb_end):
                    rfb_start = rfb_row.get("ref_start", np.nan)
                    rfb_end = rfb_row.get("ref_end", np.nan)
                    if not pd.isna(rfb_start):
                        rfb_start = (rfb_start - read_start) / 1000.0
                    if not pd.isna(rfb_end):
                        rfb_end = (rfb_end - read_start) / 1000.0
                else:
                    rfb_start = rfb_start / 1000.0
                    rfb_end = rfb_end / 1000.0

                if not pd.isna(rfb_start) and not pd.isna(rfb_end):
                    if rfb_end < rfb_start:
                        rfb_start, rfb_end = rfb_end, rfb_start
                    ax.axvline(rfb_start, color="red", linestyle="--", linewidth=1.5)
                    ax.axvline(rfb_end, color="red", linestyle="--", linewidth=1.5)

        # Stair style plot
        ax.stairs(y_vals, edges, linewidth=2, color="black", fill=False)

        # 50% BrdU probability reference line
        ax.axhline(
            y=0.50,
            color="red",
            linestyle="--",
            linewidth=1.5,
            alpha=0.8
)

        ax.set_ylim(0, 1)
        ax.set_xlabel("Position within read (kb)")
        ax.set_ylabel("BrdU probability (0–1)")
        ax.set_title(f"S_Phase Rain Plot - Read {i}\n Chromosome: {chr_label}")

        outpath = os.path.join(outdir, f"rainplot_{i:03d}_read_{i}.png")
        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)


# Runs the function that creates the rain plots
if __name__ == "__main__":
    plot_rainplots_per_read()
