# BrdU Analysis

Python script to analyze BrdU modifications in mod.bam files.

Parses through reads in a bam file from DNAscent detect .bam output and converts to a bed6 format for downstream data analysis. Uses sliding T windows. Counts number of T's in each window with a probability >= the probability given. Output is bed6 format. 

## Analyze_BrdU.py Usage
```bash
python analyze_brdu.py --bam file.mod.bam --prob 0.5 --mod ['b' 'brdu' 'e' 'edu'] --window 10 --format bed6 --output results.bed
```

## Requirements
- pysam
- numpy
