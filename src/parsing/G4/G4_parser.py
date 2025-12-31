import os
from pathlib import Path
from dotenv import load_dotenv

# Load .env to pick up the G4 extraction file path for this environment
REPO_ROOT = Path(__file__).resolve().parents[3]
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")

def parse_g4_motifs(filename):
    """
    Parse a G4 motifs text file and extract genomic coordinates.
    
    G4 motifs (G-quadruplexes) are secondary DNA structures formed by guanine-rich
    sequences. This function reads a specially formatted text file containing G4
    motif locations across chromosomes and organizes them into structured data.
    
    Expected file format:
        >chromosome_name
        Start    End    Sequence    Length
        123      456    GGGTTAGGG   10
        ...
    
    Args:
        filename (str): Path to the G4 motifs text file
        
    Returns:
        tuple: A tuple containing two dictionaries:
            - coordinates (dict): Maps chromosome names to lists of (start, end) tuples
            - lengths (dict): Maps chromosome names to lists of motif lengths
            
    Example:
        coordinates = {'chrI': [(100, 120), (300, 315)]}
        lengths = {'chrI': [20, 15]}
    """
    # Initialize empty dictionaries to store our parsed data
    coordinates = {}
    lengths = {}
    current_chr = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            # Check if the line is a chromosome header
            if line.startswith('>'):
                current_chr = line[1:] # We will remove the '>' for cleaner output
                # Initialize empty lists for this chromosome's data
                coordinates[current_chr] = []
                lengths[current_chr] = []

            # Skip empty and header lines 
            elif not line or line.startswith('Start'):
                continue

            # Parse numerical lines (Mostly (start, end) coordinates and length)
            elif current_chr is not None:
                # Split the tab-delimited line into components
                parts = line.split('\t')
                # Ensuriing we have at least four columns (Start, End, Sequence, Length)
                if len(parts) >= 4:
                    try:
                        start = int(parts[0].strip())
                        end = int(parts[1].strip())
                        motif_length = int(parts[3].strip())
                        # Storing the coordinates as a tuple
                        coordinates[current_chr].append((start, end))
                        # Store the corresponing length 
                        lengths[current_chr].append(motif_length)
                    except ValueError:
                        # Skip lines that does not have valid integers (e.g. Sequence)
                        continue

    return coordinates, lengths

# def print_chromosome_summary(coordinates, lengths):
#     '''
#     This is a test method. We will print the motifs for each chromosome
#     to the console to see if we get the expected output. 
#     '''
#     for chr_name in coordinates.keys():
#         print(f"Chromosome: {chr_name}")
#         print(f"Number of motifs: {len(coordinates[chr_name])}")
#         print(f"Coordinates: {coordinates[chr_name]}")
#         print(f"Lengths: {lengths[chr_name]}")
#         print("-" * 80)
#         print()

def save_to_bed(coordinates, lenghts, output_file):
    """
    Save G4 motif data to BED format for genome visualization tools.
    
    BED (Browser Extensible Data) is a standard tab-delimited format used by
    genome browsers and analysis tools. Each line represents a genomic feature
    with its chromosome, start, end, name, score, and strand.
    
    BED format columns:
        1. chrom: Chromosome name
        2. chromStart: Start position (0-based)
        3. chromEnd: End position (exclusive)
        4. name: Feature name
        5. score: Numeric score (here we use motif length)
        6. strand: '+', '-', or '.' for unknown
    
    Args:
        coordinates (dict): Maps chromosome names to lists of (start, end) tuples
        lengths (dict): Maps chromosome names to lists of motif lengths
        output_file (str): Path where the BED file will be written
        
    Returns:
        None: Writes output to file and prints confirmation message
        
    Example output line:
        chrI    100    120    G4_1    20    .
    """
    with open(output_file, 'w') as f:
        motif_count = 0 #Counter for creating unique motif names

        # Sort chromosomes for consistent output
        for chr_name in sorted(coordinates.keys()):
            motifs = coordinates[chr_name]
            lens = lengths[chr_name]

            # Iterate through each motif on this chromosome
            for i, (start, end) in enumerate(motifs):
                motif_count += 1
                name = f"G4_{motif_count}"
                length = lens[i]
                f.write(f"{chr_name}\t{start}\t{end}\t{name}\t{length}\t.\n")
    print(f"Saved {motif_count} motifs to {output_file}")

if __name__ == "__main__":
    g4_path = os.getenv("G4_EXTRACTION")
    if not g4_path:
        raise ValueError("G4_EXTRACTION must be set in env/.env")

    coordinates, lengths = parse_g4_motifs(g4_path)
    # print_chromosome_summary(coordinates, lengths)

    # Save to BED format, best format for plotting on bedgraph
    save_to_bed(coordinates, lengths, "g4.motifs.bed")
