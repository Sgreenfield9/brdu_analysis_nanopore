# Transposable Element (TE) Extractor (GFF to BED)

## Overview
This tool extracts transposable element (TE) genomic coordinates from GFF3 annotation files and converts them to BED format for visualization and genomic analysis.

## Features
- Parses hierarchical GFF3 annotation files with memory-efficient iteration
- Identifies TEs using multiple detection rules (feature types and qualifiers)
- Flexible output: extract TE bodies + LTRs or TE bodies only
- Handles multiple TE naming conventions
- Preserves strand and feature type information
- Extended BED format with 7 columns

## Setup
1. Set the GFF file path in `env/.env`:
```
W303_GFF=/path/to/your/annotation.gff3
```

## Input Format
Expects a standard GFF3 annotation file with TE annotations:
```
chrI    SGD    mobile_genetic_element    5000    9500    .    +    .    ID=YARCTy1-1;mobile_element_type=retrotransposon
chrI    SGD    long_terminal_repeat      5000    5330    .    +    .    ID=YARCTy1-1_LTR1;Parent=YARCTy1-1
```

## Usage
Uncomment desired option in the script, then run:
```bash
python TE_parser.py
```

**Option A - TE Bodies + LTRs (complete TE landscape):**
```python
te_bed_all = "w303_te_and_ltrs.bed"
n_all = parse_te_gff_to_bed(gff_file, te_bed_all, include_ltrs=True)
```

**Option B - TE Bodies Only (whole elements without subparts):**
```python
te_bed_bodies_only = "w303_te_bodies_only.bed"
n_bodies = parse_te_gff_to_bed(gff_file, te_bed_bodies_only, include_ltrs=False)
```

## Output
BED-formatted file with **7 columns**:
- Chromosome name
- Start position (0-based)
- End position (exclusive)
- TE name (from mobile_element_type/Name/ID)
- Feature length
- Strand (+, -, or .)
- Feature type (mobile_genetic_element or long_terminal_repeat)

## BED Format Example
```
chrI    5000    9500    retrotransposon:Ty1    4500    +    mobile_genetic_element
chrI    5000    5330    YARCTy1-1_LTR1         330     +    long_terminal_repeat
chrII   15000   19200   retrotransposon:Ty3    4200    -    mobile_genetic_element
```

## Detection Rules
**Primary Rule (Feature Type):**
- `mobile_genetic_element` → TE body
- `long_terminal_repeat` → LTR (only if `include_ltrs=True`)

**Secondary Rule (Qualifiers):**
- `gbkey="mobile_element"` → TE body
- `gbkey="repeat_region"` + `rpt_type="long_terminal_repeat"` → LTR
- Presence of `mobile_element_type` → TE body

## Notes
- Memory-efficient iteration using generators for large GFF files
- Naming priority: `mobile_element_type` → `Name` → `ID` → "TE"
- Extended BED format (7 columns) includes feature type for downstream filtering
- Works with W303 yeast strain annotations