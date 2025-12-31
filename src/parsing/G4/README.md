# G4 Motif Parser and BED Converter

## Overview
This tool parses G4-quadruplex motif data from yeast (W303 strain) and converts it into BED format for genome visualization and analysis.

## Features
- Parses G4 motif coordinates from tab-delimited text files
- Organizes data by chromosome (supports 16 chromosomes in W303)
- Exports to BED format for use with genome browsers (IGV, UCSC) and analysis tools (bedtools)
- Validates data integrity with console summary output

## Setup
1. Set the input file path in `env/.env`:
```
   G4_EXTRACTION=/path/to/your/g4_motifs.txt
```

## Input Format
Expected input file structure:
```
>chrI
Start    End    Sequence    Length
123      145    GGGTTAGGG   22
...
>chrII
Start    End    Sequence    Length
...
```

## Usage
```bash
python G4_parser.py
```

## Output
- **Console**: Summary of motifs per chromosome
- **g4.motifs.bed**: BED-formatted file with columns:
  - Chromosome name
  - Start position
  - End position
  - Motif name (G4_1, G4_2, ...)
  - Length (score field)
  - Strand (.)

## BED Format Example
```
chrI    123    145    G4_1    22    .
chrI    300    318    G4_2    18    .
chrII   450    472    G4_3    22    .
```

## Use Case
- Find overlaps with other genomic features using **bedtools**
