# BrdU Analysis

Python script to analyze BrdU modifications in mod.bam files.

## Usage
```bash
python analyze_brdu.py --bam file.mod.bam --prob 0.5 --mod ['b' 'brdu' 'e' 'edu'] --window 10 --format bed6 --output results.bed
```

## Requirements
- pysam
- numpy
