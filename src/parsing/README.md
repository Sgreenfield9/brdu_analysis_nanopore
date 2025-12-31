# Parsing

## About
We have three different parsers. We will have a parser for G4 motifs, tRNA, and Transposable Elements. The G4 parser extracts start/end coordinates along with the length from a .txt file. We got the .txt file from using the tool G4Hunter. Both the tRNA and Transposable elements will using BioPython GFF to parse the GFF files. We obtained the GFF files from NBCI using there datasets command-line-tools. 

## How to run the parsers

### G4 Parser
If you are running from the root directory:
```bash
python3 src/parsing/G4/G4_parser.py
```

### tRNA Parser
If you are running from the root directory:
```bash
python3 src/parsing/tRNA/tRNA_parser.py
```

### Transposable Elements
Still working on this method 

## Details
More details about how these parsers work will be provided in a markdown under their designated directory. We will have a G4, TE (Transposable Elements), and tRNA directory each containing their parser. 


