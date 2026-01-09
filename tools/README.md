This is the tools page. This is where I describe how to use the open-source tools for this project.

# G4Hunter

G4Hunter is a tool that helps extract G4's from .fasta file formats. Thanks to [AnimaTardeb](https://github.com/AnimaTardeb/G4Hunter) for making this tool. These are the steps you need to run in order to replicate my steps.

1. Create the Python 2.7 virtual environment in the **root directory**
    ```bash
    conda create \
  --prefix ./g4hunter_env \
  python=2.7 \
  biopython \
  numpy
  ```
2. Then you will activate the virtual environment
```bash
    conda activate ./env/g4hunter_env/
```
3. Go to the directory that contains G4Hunter
```bash
    cd tools/G4Hunter_env
```
4. Run the following command to get and output 
```bash
    python G4Hunter.py   -i ../../data/ncbi/ncbi_dataset/GCF_000146045.2_R64_genomic.fna   -o ../../results/G4_extraction   -w 25   -s 1.5
```
5. Deactivate this virtual environment
```bash
    conda deactivate
```

# TelFinder

TelFinder is a tool that helps find subtelomeric regions. But, the tool is more used to find scores and give details on what makes it a subtelomeric region. This tool is not necessary for this project, but if we ever want to use it examples of how to run it will be provided below.
```bash
    python3 TelFinder.py -f fa -inf /home/kylep/MAPs/prj_01/brdu_analysis_nanopore/data/ncbi/ncbi_dataset/GCF_000146045.2_R64_genomic.fna -s NCR -o /home/kylep/MAPs/prj_01/brdu_analysis_nanopore/results/subtelomeric/left -e left
```
```bash
    python3 TelFinder.py -f fa -inf /home/kylep/MAPs/prj_01/brdu_analysis_nanopore/data/ncbi/ncbi_dataset/GCF_000146045.2_R64_genomic.fna -s NCR -o /home/kylep/MAPs/prj_01/brdu_analysis_nanopore/results/subtelomeric/right -e right
```

The GitHUb repository that has the source code for this tool can be found [here](https://github.com/bio-tombow/TelFinder).