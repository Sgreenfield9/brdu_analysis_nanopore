# BrdU Percentage Exporter

## Overview
This module computes BrdU incorporation percentage from combined bedgraph data and exports sorted CSVs for downstream plotting and analysis.

## Features
- Calculates BrdU% as `(BrdU_count / Nmod) * 100` with smoothing
- Deduplicates positions per chromosome and optionally enforces minimum spacing
- Exports sorted CSVs (highest BrdU% first) for M phase and S phase

## Setup
1. Set output paths in `env/.env`:
```
BRDU_PCT_OUTPUT_M=/path/to/output/BrdU_pct/M_Phase
BRDU_PCT_OUTPUT_S=/path/to/output/BrdU_pct/S_Phase
```

## Usage
```bash
python src/parsing/BrdU_pct/BrdU_pct_M.py
python src/parsing/BrdU_pct/BrdU_pct_S.py
```

## Output
- **BrdU_pct_sorted_desc_2.5kb.csv**
- **BrdU_pct_sorted_desc_5kb.csv**
- **BrdU_pct_sorted_desc_10kb.csv**

Each CSV includes: `chrom, start, end, BrdU_pct, BrdU_count, Nmod`.