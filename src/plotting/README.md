# Plotting

## Overview
Plotting scripts for BrdU incorporation across yeast chromosomes in M phase and S phase. Outputs include genome browser tracks, point of interest windows, rain plots from single reads, and probability threshold diagnostics.

### M Phase Plots
#### Genome Browsers
Full chromosome browser plots for M phase with two variants. The smoothed version uses rolling windows to display Nmod, BrdU counts, and BrdU percent across entire chromosomes. The unsmoothed version uses raw signals. Both place genomic feature rulers below the main track, including G4, tRNA, TE, and subtelomeric regions.
Scripts: 
- `M_phase_chromosome_plotting.py`
- `M_phase_chromosome_plotting_unsmoothed.py`.

#### POI Plots (Top 20 from Genome Browsers)
Point of interest windows for the top 20 BrdU percent sites identified from the genome browser ranking files. Each plot centers a 30 kb window on the site and shows the same track layout as the genome browser. Smoothed and unsmoothed variants are provided to compare denoised versus raw signals.
Scripts:
- `M_phase_chromosome_plotting_pct_10kb.py`
- `M_phase_chromosome_plotting_pct_unsmoothed.py`.

#### Rain Plots
Single read rain plots for M phase that estimate BrdU incorporation probability using a sliding window over T bases. The scatter points are per base BrdU probability and the stair line is the mean within each T window. There is a smoothed version and a no smoothing version for comparison.
Scripts: 
- `rain_plots_single_read.py`
- `rain_plots_single_read_unsmoothed.py`.

### S Phase Plots
#### Genome Browsers
Full chromosome browser plots for S phase with smoothed and unsmoothed variants. Tracks include Nmod, BrdU counts, BrdU percent, and feature rulers for G4, tRNA, TE, and subtelomeric regions.
Scripts:
- `S_phase_chromosome_plotting.py`
- `S_phase_chromosome_plotting_unsmoothed.py`.

#### POI Plots (Top 20 from Genome Browsers)
Point of interest windows for the top 20 BrdU percent sites identified from the genome browser ranking files. Each plot centers a 30 kb window on the site and shows the same track layout as the genome browser. Smoothed and unsmoothed variants are provided to compare denoised versus raw signals.
Scripts: 
- `S_phase_chromosome_plotting_pct.py`
- `S_phase_chromosome_plotting_pct_unsmoothed.py`.

#### Rain Plots
Single read and rDNA focused rain plots for S phase that estimate BrdU incorporation probability using T based sliding windows. The scatter points are per base BrdU probability and the stair line is the mean within each T window. rDNA variants add RFB overlays and sliding window thresholding, plus a multi read averaged stair plot.
Scripts:
- `rain_plots_single_read.py`
- `rain_plots_single_read_unsmoothed.py`
- `rain_plots_single_read_rDNA.py`, 

### Probability Threshold Plots
#### Probability Thresholds
Diagnostics for ML probability cutoffs. One script compares positive and negative distributions across thresholds, the other shows positive control boxplots including unfiltered scores.
Scripts: 
- `plot_ml_distributions.py`
- `plott_thresholds.py`.
