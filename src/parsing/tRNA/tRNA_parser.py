"""
tRNA Feature Extractor from GFF to BED Format

This script extracts tRNA (transfer RNA) features from a GFF3 annotation file
and converts them into BED format for genome visualization and analysis.

tRNAs are essential RNA molecules that decode mRNA sequences into amino acids
during protein synthesis. This tool helps visualize their genomic locations.
"""

import os
from BCBio import GFF
from dotenv import dotenv_values

def iter_features(features):
    """
    Recursively iterate through all features and their sub-features in a GFF record.
    
    GFF files have a hierarchical structure where features can contain sub-features
    (e.g., a gene can contain mRNA, which contains exons). This generator function
    flattens this hierarchy to allow iteration over all features at all levels.
    
    Args:
        features (list): List of SeqFeature objects from a GFF record
        
    Yields:
        SeqFeature: Each feature and all of its nested sub-features
        
    Example:
        Gene
        ├── mRNA
        │   ├── exon
        │   └── CDS
        └── tRNA
        
        This function yields: Gene, mRNA, exon, CDS, tRNA
    """
    for feature in features:
        yield feature
        # Check if this feature has an sub-features
        sub_features = getattr(feature, "sub_features", None)
        # If sub-feature exists, recursively iterate through them
        if sub_features:
            for sub_feature in iter_features(sub_features):
                yield sub_feature

# Construct the path to the root directory
# Go up three levels (".." x3) to get to script location
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
# Path to the .env file containing environment variables
env_path = os.path.join(repo_root, "env", ".env")
# Load environment variables from the .env file
env_vars = dotenv_values(env_path)
# Retrieve the path to the W303 GFF annotation file
gff_file = env_vars.get("W303_GFF")
if not gff_file:
    raise ValueError("W303_GFF not found in env/.env")

bed_file = "trna_coordinates.bed"

with open(gff_file) as in_handle, open(bed_file, "w") as out_handle:
    """
    Extract tRNA features from a GFF file and write them to BED format.
    
    This function reads a GFF3 annotation file, filters for tRNA features,
    and converts their genomic coordinates and metadata into BED format
    suitable for genome browsers and analysis tools.
    
    Process:
        1. Parse the GFF file to access all chromosomes/contigs
        2. Recursively iterate through all features (including nested ones)
        3. Filter for tRNA-type features only
        4. Extract: chromosome, coordinates, name, and strand information
        5. Write each tRNA as a BED-formatted line
    
    Args:
        gff_file (str): Path to the input GFF3 annotation file
        bed_file (str): Path to the output BED file to create
        
    Returns:
        None: Writes output directly to the BED file
        
    BED Output Format:
        Column 1: Chromosome name
        Column 2: Start position (0-based)
        Column 3: End position (exclusive)
        Column 4: tRNA name (gene/Name/ID)
        Column 5: Feature length
        Column 6: Strand (+, -, or .)
    """
    for record in GFF.parse(in_handle):
        for feature in iter_features(record.features):
            if feature.type != "tRNA":
                continue

            chrom = record.id
            start = int(feature.location.start)
            end = int(feature.location.end)
            length = end - start

            name = feature.qualifiers.get("gene", [None])[0]
            if not name:
                name = feature.qualifiers.get("Name", [None])[0]
            if not name:
                name = feature.qualifiers.get("ID", ["tRNA"])[0]

            strand_val = feature.location.strand
            strand = "+" if strand_val == 1 else "-" if strand_val == -1 else "."

            out_handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\n")

print(f"BED file created: {bed_file}")
