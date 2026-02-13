import os
from pathlib import Path
from dotenv import load_dotenv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[4]
load_dotenv(dotenv_path=REPO_ROOT /"env" / ".env")

def load_rain_plot_input():
    """
    Finds our rainplot input file path. It will then convert
    that file into a pandas dataframe to make it easier/faster to plot.
    """
    rain_plot_path = os.getenv('INPUT_FILE_RP_S')

    if not rain_plot_path:
        raise ValueError("INPUT_FILE_RP must be set in the env/.env")
    
    rain_plot_df = pd.read_csv(
        rain_plot_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "read_id", "mod_base", "mod_prob"]
    )

    return rain_plot_df

def get_output_dir():
    """
    Determine the output directory for the rain plots
    """
    output_dir = os.getenv("OUTPUT_FILE_RP_UNSMOOTHED_S")

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
    Plot rain plots for 100 reads with the raw data:

    We use the raw data, meaning there is no smoothing.
    Iterate through all the reads, assign them a group id,
    assign them a chromosome id, and plot the graph with scatter 
    points and stairs line.
    """
    df = load_rain_plot_input()
    outdir = get_output_dir()

    # small memory optimization
    df["start"] = pd.to_numeric(df["start"], errors="coerce", downcast="integer")
    df["mod_prob"] = pd.to_numeric(df["mod_prob"], errors="coerce", downcast="float")
    df = df.dropna(subset=["read_id", "start", "mod_prob"])

    max_reads = 100

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

        read_start = sub["start"].min()
        x = (sub["start"].to_numpy() - read_start) / 1000.0  # kb
        y = sub["mod_prob"].to_numpy()

        # Get the number of bases in a genomic read
        # If less than two we will skip 
        if len(x) < 2:
            continue

        # Build edges without smoothing
        # edges are midpoints between consecutive x's, with an extra edge on each side.
        # This preserves every raw y value as a step height.
        x = x.astype(float)
        y = y.astype(float)

        mids = (x[:-1] + x[1:]) / 2.0
        left_edge = x[0] - (mids[0] - x[0])
        right_edge = x[-1] + (x[-1] - mids[-1])

        edges = np.concatenate(([left_edge], mids, [right_edge]))

        # Ensure edges are non-decreasing (just in case of duplicate positions)
        edges = np.maximum.accumulate(edges)

        fig, ax = plt.subplots(figsize=(10, 4))

        # Raw rain scatter (keep jumbled)
        ax.scatter(x, y, s=6, alpha=0.5, color="black", zorder=1)

        # Drawn as steps (will look the same as line because there is no smoothing)
        ax.stairs(y, edges, linewidth=1, alpha=0.9, color="black", fill=False, zorder=2)

        ax.set_ylim(0, 1)
        ax.set_xlabel("Position within read (kb)")
        ax.set_ylabel("BrdU probability (0–1)")
        ax.set_title(f"S_Phase (Raw Data) Rain Plot - Read {i} \n Chromosome: {chr_label}")
        ax.grid(True, linewidth=0.3, alpha=0.4)

        outpath = os.path.join(outdir, f"rainplot_{i:03d}_read_{i}_unsmoothed.png")

        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)


# Runs the rainplot function to generate plots
if __name__ == "__main__":
    plot_rainplots_per_read()
