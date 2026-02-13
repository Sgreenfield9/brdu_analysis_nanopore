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
    rain_plot_path = os.getenv("INPUT_FILE_RP")

    if not rain_plot_path:
        raise ValueError("INPUT_FILE_RP must be set in the env/.env")

    rain_plot_df = pd.read_csv(
        rain_plot_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "read_id", "mod_base", "mod_prob"],
    )

    return rain_plot_df


def get_output_dir():
    """
    Determine the output directory for the rain plots
    """
    output_dir = os.getenv("OUTPUT_FILE_RP")

    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "output", "rain_plot")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

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

    # Binning control
    points_per_bin = 200  # Can change to 100 for more detail, 500 for smoother plots
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


        # Stair style plot
        ax.stairs(y_vals, edges, linewidth=2, color="black", fill=False)

        ax.set_ylim(0, 1)
        ax.set_xlabel("Position within read (kb)")
        ax.set_ylabel("BrdU probability (0–1)")
        ax.set_title(f"Mitosis Rain Plot - Read {i}")

        outpath = os.path.join(outdir, f"rainplot_{i:03d}_read_{i}.png")
        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)


# Runs the function that creates the rain plots
if __name__ == "__main__":
    plot_rainplots_per_read()
