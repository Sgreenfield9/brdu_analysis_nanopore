# tRNA Feature Extractor (GFF to BED)

## Overview
This tool extracts tRNA (transfer RNA) genomic coordinates from GFF3 annotation files and converts them to BED format for visualization and analysis.

## Features
- Parses hierarchical GFF3 annotation files
- Recursively extracts tRNA features from nested annotations
- Handles multiple tRNA naming conventions (gene, Name, ID)
- Preserves strand information (+/-)
- Exports to BED format for genome browser compatibility

## Setup
1. Set the GFF file path in `env/.env`:
```
   W303_GFF=/path/to/your/annotation.gff3
```

## Input Format
Expects a standard GFF3 annotation file with tRNA features:
```
chrI    SGD    tRNA    1000    1072    .    +    .    ID=tRNA:tA(AGC)A;Name=tA(AGC)A;gene=tA(AGC)A
chrI    SGD    tRNA    2500    2571    .    -    .    ID=tRNA:tD(GUC)B;Name=tD(GUC)B;gene=tD(GUC)B
```

## Usage
```bash
python tRNA_parser.py
```

## Output
**trna_coordinates.bed**: BED-formatted file with columns:
- Chromosome name
- Start position (0-based)
- End position (exclusive)
- tRNA name (extracted from gene/Name/ID fields)
- Feature length
- Strand (+, -, or .)

## BED Format Example
```
chrI    1000    1072    tA(AGC)A    72    +
chrI    2500    2571    tD(GUC)B    71    -
chrII   3400    3471    tE(UUC)C    71    +
```

## Use Case
- Find overlaps with other genomic features using **bedtools**

