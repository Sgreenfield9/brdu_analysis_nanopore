This is the tools page. This is where I describe how to use the open-source tools for this project.

# G4Hunter

G4Hunter is a tool that helps extract G4's from Fasta file formats. Thanks to [AnimaTardeb](https://github.com/AnimaTardeb/G4Hunter) for making this tool. These are the steps you need to run in order to replicate my steps.

1. Create the Python 2.7 virtual environment in the **root directory**
    conda create \
  --prefix ./g4hunter_env \
  python=2.7 \
  biopython \
  numpy
2. Then you will activate the virtual environment
    conda activate ./env/g4hunter_env/
3. Go to the directory that contaings G4Hunter
    cd tools/G4Hunter_env
4. Run the following command to get and output 
    python G4Hunter.py   -i ../../data/ncbi/ncbi_dataset/GCF_000146045.2_R64_genomic.fna   -o ../../results/G4_extraction   -w 25   -s 1.5
5. Deactivate this virtual environment
    conda deactivate