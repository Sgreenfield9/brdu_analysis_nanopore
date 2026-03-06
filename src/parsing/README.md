# Parsing

## About
This is the parsing directory which will contain all the parsers needed to extract data from public datasets or CU Anschutz datasets. We currently have three parsers for genomic features (G4's, tRNA, and Transposable Elements). I used a tool called G4 Hunter (more about this under the tools directory) which creates a .txt file and our parser will extract those well known G4 motifs. Both the tRNA and Transposable Elements will be using BioPython GFF (an open-source library) to parse through the publicly known NCBI .GFF dataset. Lastly, there is a parser for BrdU_pct parser for both M Phase & S Phase. This helps us identify our top 20 points of interest in our dataset. 

## How to run the parsers

### G4 Parser
If you are running from the root directory:
```bash
python3 src/parsing/G4/G4_parser.py
```
or under this directory:
```bash
python3 G4/G4_parser.py
```

### tRNA Parser
If you are running from the root directory:
```bash
python3 src/parsing/tRNA/tRNA_parser.py
```
or under this directory:
```bash
python3 tRNA/tRNA_parser.py
```

### Transposable Elements
If you are running from the the root directory:
```bash
python3 src/parsing/TE/TE_parser.py
```
or under this directory:
```bash
python3 TE/TE_parser.py
```

### BrdU_pct
If you are running from the root directory:
```bash
python3 src/parsing/BrdU_pct/BrdU_pct_M.py
python3 src/parsing/BrdU_pct/BrdU_pct_S.py
```
or under this directory:
```bash
python3 BrdU_pct/BrdU_pct_M.py
python3 BrdU_pct/BrdU_pct_S.py
```

## Details
More details about how these parsers work will be provided in a markdown under their designated directory. We will have a G4, TE (Transposable Elements), BrdU pct, and tRNA directory each containing their parser. 