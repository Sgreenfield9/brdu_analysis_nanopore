# BrdU Analysis of Nanopore 
## Requirements
First create a virtual environment by running these commands:
```bash
virtualenv .venv
source .venv/bin/activate
```
Once you have your virtual environment created and running, run the following command to install the requirements:
```bash
pip install -r requirements.txt
```
## Overview
This project provides a pipeline for extracting and visualizing BrdU related genomic features in Nanopore sequencing data for the W303 strain on yeast. It uses publicly available NCBI genome data to parse G4 motifs, tRNA, and transposable elements, then plots them on a chromosome browser.
## Data Acquisition
**NOTE: Make sure to install NCBI Datasets CLI tool. The documentation is found [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/):**
```bash
datasets download genome accession GCA_002163515.1 \
 --include gff3,genome
```
```bash
unzip ncbi_dataset.zip
``` 
## Parsing
Parsing scripts in `src/parsing` extract genome feature tracks and prepare BrdU related inputs used by the plotting pipeline.

## Plotting
Plotting scripts in `src/plotting` generate genome browser plots, point of interest windows, and rain plots for BrdU incorporation.

## Results
**Work in progress**

## References/Resources
### NCBI
This project uses publicly available genomic data from National Center for Biotechnology Information ([NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_002163515.1/)):

- Organism : *Saccharomyces cerevisiae* (W303)

Data was downloaded using the NCBI Datasets CLI, accessed December 2025.
### CU Anschutz (McClure Lab)
BrdU sequencing datasets were provided by CU Anschutz, McClure Lab.
### Biopython
This project uses the Biopython tool to extract .GFF file embeddings for important genome features. It is also used to generate the genome diagrams. The documentation can be found [here](https://biopython.org/). In this project we used version 1.86. 
### G4 Hunter Tool 
This project uses the G4 Hunter tool created by **AnimaTardeb** to extract G4 motifs from a .fasta file format. The latest version of this tool is found on [GitHub](https://github.com/AnimaTardeb/G4Hunter).
