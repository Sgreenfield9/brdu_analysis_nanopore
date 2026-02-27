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


def get_output_dir():
    """
    Determine the output directory for the rain plots
    """
    output_dir = os.getenv("OUTPUT_FILE_RP_rDNA_SW")

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
       - We now compute the stairs using a sliding window of 100 T bases.
       - For each window we calculate the proportion of T bases with mod_prob > 0.50.
       - That proportion (0 to 1) is what we plot as the "stairs".

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

    min_T_count = 1500 # change this number to whatever threshold you want

    # This is the code that makes it random each time you run the script.
    # Comment this block out if you want the first 100 reads every time.
    rng = np.random.default_rng()

    # Normalize mod_base in case of lowercase / weird formatting
    df["mod_base"] = df["mod_base"].astype(str).str.upper()


    window_T = 75            # Sliding window size in number of T bases
    prob_thresh = 0.50       # "above 50%" threshold for mod_prob

    # Instead of overlaying multiple read stairs, we compute ONE averaged stair
    # per group of reads by binning each read onto a common kb-axis and averaging.
    reads_per_plot = 5       # Change this number to control how many reads are averaged per plot
    max_plots = 50           # Change this number to generate more than 20 plots

    # Common kb bins (this controls the "stair resolution" of the averaged curve)
    kb_bin = 0.10            # 0.10 kb = 100 bp bins (increase to 0.2/0.5 for smoother)
    max_kb = 12.0            # Maximum x-range to support (in kb) for averaging across reads

    # Require that the averaged curve has at least one bin >= 0.05
    # (prevents averaged flat plots)
    min_avg_peak = 0.05

    # First: compute T counts per read, then filter to reads that pass threshold
    per_read = (
        df.groupby("read_id", as_index=False)
          .agg(
              T_count=("mod_base", lambda s: (s == "T").sum())
          )
    )

    min_T_count = max(min_T_count, window_T)

    eligible_ids = per_read.loc[per_read["T_count"] >= min_T_count, "read_id"].to_numpy()

    if eligible_ids.size == 0:
        raise ValueError(
            f"No reads passed the T threshold. min_T_count={min_T_count}. "
            f"Lower min_T_count or verify mod_base contains 'T'."
        )

    # Then: filter those eligible reads so we do not get flat lines
    # We keep your variability checks but now they operate on the per-read rolling proportion signal.
    min_prop_spread = 0.10
    min_prop_range = 0.10
    min_prop_peak = 0.05

    df_eligible = df[df["read_id"].isin(eligible_ids)]
    keep_ids = []

    for rid, sub in df_eligible.groupby("read_id", sort=False):
        sub = sub.sort_values("start", kind="mergesort")

        sub_T = sub[sub["mod_base"] == "T"]
        if len(sub_T) < window_T:
            continue

        yT = sub_T["mod_prob"].to_numpy(dtype=float)
        above = (yT > prob_thresh).astype(np.int32)

        csum = np.cumsum(above, dtype=np.int64)
        window_sum = csum[window_T - 1:] - np.concatenate(([0], csum[:-window_T]))
        prop = window_sum / float(window_T)

        if prop.size < 2:
            continue

        # Ensure each read has at least one window with some signal
        if float(prop.max()) < min_prop_peak:
            continue

        q05 = float(np.quantile(prop, 0.05))
        q95 = float(np.quantile(prop, 0.95))
        prop_spread = q95 - q05
        if prop_spread < min_prop_spread:
            continue

        prop_range = float(prop.max() - prop.min())
        if prop_range < min_prop_range:
            continue

        keep_ids.append(rid)

    keep_ids = np.asarray(keep_ids, dtype=object)

    if keep_ids.size == 0:
        raise ValueError(
            f"No reads passed the variability filter. "
            f"Try lowering min_prop_spread={min_prop_spread} or min_prop_range={min_prop_range} "
            f"or min_prop_peak={min_prop_peak}."
        )

    # Sample up to max_reads from those eligible reads
    if keep_ids.size > max_reads:
        sampled_ids = rng.choice(keep_ids, size=max_reads, replace=False)
    else:
        sampled_ids = keep_ids

    # Keep only sampled reads
    df = df[df["read_id"].isin(sampled_ids)]

    # Optional scatter: we will plot the raw points from ALL reads in the group (like your old plots)
    show_scatter = True

    # Create common bin edges for averaging (stairs will use these edges)
    # Example: kb_bin=0.1 -> edges at 0.0, 0.1, 0.2, ... max_kb
    bin_edges = np.arange(0.0, max_kb + kb_bin, kb_bin, dtype=float)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    # We will walk reads in order and process in groups of reads_per_plot.
    grouped = list(df.groupby("read_id", sort=False))

    plot_count = 0
    group_start_index = 0

    while group_start_index < len(grouped) and plot_count < max_plots:
        # Take the next chunk of reads
        group_slice = grouped[group_start_index: group_start_index + reads_per_plot]
        if len(group_slice) == 0:
            break

        # NEW CODE:
        # We will collect per-read binned "prop_above" arrays and then average them.
        per_read_binned = []
        all_scatter_x = []
        all_scatter_y = []
        chr_labels = []

        for local_i, (_, sub) in enumerate(group_slice, start=1):
            sub = sub.sort_values("start", kind="mergesort")

            genbank_ids = sub["chrom"].astype(str).unique()
            if len(genbank_ids) == 1:
                chr_label = genbank_to_chr.get(genbank_ids[0], genbank_ids[0])
            else:
                chr_label = "mixed/" + ",".join(genbank_to_chr.get(g, g) for g in genbank_ids)
            chr_labels.append(str(chr_label))

            read_start = sub["start"].min()
            x = (sub["start"].to_numpy() - read_start) / 1000.0
            y = sub["mod_prob"].to_numpy(dtype=float)

            # Collect scatter points from this read (optional)
            if show_scatter:
                all_scatter_x.append(x)
                all_scatter_y.append(y)

            sub_T = sub[sub["mod_base"] == "T"]
            if len(sub_T) < window_T:
                continue

            xT = (sub_T["start"].to_numpy() - read_start) / 1000.0
            yT = sub_T["mod_prob"].to_numpy(dtype=float)

            above = (yT > prob_thresh).astype(np.int32)
            csum = np.cumsum(above, dtype=np.int64)
            window_sum = csum[window_T - 1:] - np.concatenate(([0], csum[:-window_T]))
            prop_above = window_sum / float(window_T)

            # Each rolling window corresponds to a window start index in xT.
            # We map each window's prop_above to the kb position of its window-start,
            # then bin onto our common kb bins so different reads can be averaged.
            window_start_x = xT[:prop_above.size]

            # Bin the window_start_x values into fixed kb bins
            # We compute mean prop_above within each bin.
            bin_idx = np.digitize(window_start_x, bin_edges) - 1  # bins 0..len-2
            valid = (bin_idx >= 0) & (bin_idx < (len(bin_edges) - 1))
            if not np.any(valid):
                continue

            tmp = pd.DataFrame({"bin": bin_idx[valid], "p": prop_above[valid]})
            binned = tmp.groupby("bin", as_index=True)["p"].mean()

            # Build a full-length array (NaN where this read has no data in that bin)
            y_bins = np.full(len(bin_edges) - 1, np.nan, dtype=float)
            y_bins[binned.index.to_numpy(dtype=int)] = binned.to_numpy(dtype=float)

            per_read_binned.append(y_bins)

        # If we didn't get enough reads with valid binned curves, skip this group
        if len(per_read_binned) == 0:
            group_start_index += reads_per_plot
            continue

        # Average across reads per bin, ignoring NaNs
        stacked = np.vstack(per_read_binned)
        avg_curve = np.nanmean(stacked, axis=0)

        # If the averaged curve is basically flat (no signal), skip saving this plot
        # (also helps avoid the "all zeros" look)
        if not np.isfinite(avg_curve).any() or float(np.nanmax(avg_curve)) < min_avg_peak:
            group_start_index += reads_per_plot
            continue

        # Convert avg_curve into stairs values/edges by dropping bins that are all NaN
        valid_bins = np.isfinite(avg_curve)
        if valid_bins.sum() < 2:
            group_start_index += reads_per_plot
            continue

        # We will keep only contiguous bins that have data
        # (If you prefer to keep gaps, we can fill NaNs with 0 instead.)
        kept_bins = np.where(valid_bins)[0]
        first_bin = int(kept_bins[0])
        last_bin = int(kept_bins[-1])

        y_vals = avg_curve[first_bin:last_bin + 1]
        edges = bin_edges[first_bin:last_bin + 2]  # +2 because edges are one longer than values

        # Ensure edges are non-decreasing
        edges = np.maximum.accumulate(edges)

        fig, ax = plt.subplots(figsize=(16, 4))

        # Optional scatter: plot all raw points from the group
        if show_scatter:
            xs = np.concatenate(all_scatter_x) if len(all_scatter_x) else np.array([], dtype=float)
            ys = np.concatenate(all_scatter_y) if len(all_scatter_y) else np.array([], dtype=float)
            if xs.size and ys.size:
                ax.scatter(xs, ys, s=1, color="black", alpha=0.15)

        # Single averaged stair
        ax.stairs(y_vals, edges, linewidth=2, color="black", fill=False)

        ax.set_ylim(0, 1)
        ax.set_xlabel("Position within read (kb)")
        ax.set_ylabel(f"Avg proportion of T's with BrdU prob > {prob_thresh:.2f}\n(rolling window = {window_T} T bases)")
        
        uniq_chr = []
        for c in chr_labels:
            if c not in uniq_chr:
                uniq_chr.append(c)
        chr_summary = ", ".join(uniq_chr[:4])
        if len(uniq_chr) > 4:
            chr_summary = chr_summary + ", ..."

        start_read_num = group_start_index + 1
        end_read_num = group_start_index + len(group_slice)

        ax.set_title(
            f"S_Phase Rain Plot - Reads {start_read_num} to {end_read_num}\n"
            f"Chromosome(s): {chr_summary} (SW-AVG)"
        )

        plot_count += 1
        outpath = os.path.join(outdir, f"rainplot_avg_{plot_count:03d}_reads_{start_read_num:03d}_{end_read_num:03d}.png")
        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)

        group_start_index += reads_per_plot


# Runs the function that creates the rain plots
if __name__ == "__main__":
    plot_rainplots_per_read()