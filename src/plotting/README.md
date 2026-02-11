# BrdU Plotting 

## Overview
This directory contains plotting scripts that visualize BrdU incorporation across the yeast genome for M phase and S phase. The plots combine smoothed long read coverage, BrdU counts, BrdU percentage, and genomic feature tracks (G4 motifs, tRNA genes, and transposable elements), with subtelomeric regions highlighted for context.

## Features
- Generates per-chromosome genome browser plots for M phase and S phase
- Computes and visualizes smoothed Coverage (Nmod), BrdU counts, and BrdU% tracks
- Overlays genomic features (G4, tRNA, TE) plus subtelomeric region tracks
- Produces focused 10 kb windows centered on top BrdU% locations from precomputed CSVs

## Setup
1. Set required paths in `env/.env` (examples shown below):
```
POSITIVE_BEDGRAPH_M=/path/to/m_phase_positive.bedgraph
NEGATIVE_BEDGRAPH_M=/path/to/m_phase_negative.bedgraph
POSITIVE_BEDGRAPH_S=/path/to/s_phase_positive.bedgraph
NEGATIVE_BEDGRAPH_S=/path/to/s_phase_negative.bedgraph
G4_MOTIFS_BED=/path/to/g4.motifs.bed
tRNA_MOTIFS_BED=/path/to/trna.motifs.bed
TE_MOTIFS_BED_ALL=/path/to/te.motifs.bed
OUTPUT_DIR_M=/path/to/output/m_phase
OUTPUT_DIR_S=/path/to/output/s_phase
BRDU_M_PCT_10KB=/path/to/BrdU_pct_sorted_desc_10kb.csv
BRDU_M_PCT_10KB_OUTPUT=/path/to/output/m_phase_10kb
BRDU_S_PCT_10KB=/path/to/BrdU_pct_sorted_desc_10kb.csv
BRDU_S_PCT_10KB_OUTPUT=/path/to/output/s_phase_10kb
```

## Scripts
- **genome_browser/M_phase_chromosome_plotting.py**
  - Produces full-chromosome plots for M phase, including smoothed Coverage, BrdU counts, BrdU%, feature tracks, and subtelomeric regions.
- **genome_browser/S_phase_chromosome_plotting.py**
  - Produces full-chromosome plots for S phase with the same tracks and styling as M phase.
- **POI_plots/M_phase_chromosome_plotting_pct_10kb.py**
  - Generates 10 kb window plots centered on the top 20 BrdU% locations from the M-phase CSV.
  - The BrdU% ruler uses the CSV value for the window center.
- **POI_plots/S_phase_chromosome_plotting_pct.py**
  - Generates 10 kb window plots centered on the top 20 BrdU% locations from the S-phase CSV.
  - The BrdU% ruler uses the CSV value for the window center.

## Usage
```bash
python src/plotting/M_Phase/genome_browser/M_phase_chromosome_plotting.py
python src/plotting/S_Phase/genome_browser/S_phase_chromosome_plotting.py
python src/plotting/M_Phase/POI_plots/M_phase_chromosome_plotting_pct_10kb.py
python src/plotting/S_Phase/POI_plots/S_phase_chromosome_plotting_pct.py
```

## Output
- **Full-chromosome plots**: Saved to `OUTPUT_DIR_M` and `OUTPUT_DIR_S`.
- **10 kb window plots**: Saved to `BRDU_M_PCT_10KB_OUTPUT` and `BRDU_S_PCT_10KB_OUTPUT`.
- Each PNG includes:
  - Main coverage/BrdU pileup track
  - Ruler panel for Coverage, BrdU, BrdU%, and genomic features
  - Subtelomeric regions highlighted in red

## Notes
- BrdU% is computed as `(BrdU_count / Nmod) * 100` and smoothed with a rolling window for visualization.
- The 10 kb window plots use the precomputed CSV ordering for top location and zoomed-in visualization.
