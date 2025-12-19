import os
from pathlib import Path
from dotenv import load_dotenv

# Load .env to pick up the G4 extraction file path for this environment
REPO_ROOT = Path(__file__).resolve().parents[3]
load_dotenv(dotenv_path=REPO_ROOT / "env" / ".env")

def parse_g4_motifs(filename):
    '''
    Parse the G4 motifs .txt file. We will store the start/end/length coordinates
    for each chromosome. There should be 16 chromomsomes in a W303 strand.

    Create two dictionaries:
     - coordinates: Will contain the chromosome start/end coordinates
     - length: Will contatin the length bewtween the start/end coordinates 
    '''
    coordinates = {}
    lengths = {}
    current_chr = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            # Check if the line is a chromosome header
            if line.startswith('>'):
                current_chr = line[1:] # We will remove the '>' for cleaner output
                coordinates[current_chr] = []
                lengths[current_chr] = []

            # Skip empty and header lines 
            elif not line or line.startswith('Start'):
                continue

            # Parse numerical lines (Mostly (start, end) coordinates and length)
            elif current_chr is not None:
                parts = line.split('\t')
                if len(parts) >= 4:
                    try:
                        start = int(parts[0].strip())
                        end = int(parts[1].strip())
                        motif_length = int(parts[3].strip())
                        coordinates[current_chr].append((start, end))
                        lengths[current_chr].append(motif_length)
                    except ValueError:
                        # Skip lines that does not have valid integers (e.g. Sequence)
                        continue

    return coordinates, lengths

def print_chromosome_summary(coordinates, lengths):
    '''
    This is a test method. We will print the motifs for each chromosome
    to the console to see if we get the expected output. 
    '''
    for chr_name in coordinates.keys():
        print(f"Chromosome: {chr_name}")
        print(f"Number of motifs: {len(coordinates[chr_name])}")
        print(f"Coordinates: {coordinates[chr_name]}")
        print(f"Lengths: {lengths[chr_name]}")
        print("-" * 80)
        print()

def save_to_bed(coordinates, lenghts, output_file):
    with open(output_file, 'w') as f:
        motif_count = 0
        for chr_name in sorted(coordinates.keys()):
            motifs = coordinates[chr_name]
            lens = lengths[chr_name]
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
    print_chromosome_summary(coordinates, lengths)

    # Save to BED format, best format for plotting on bedgraph
    save_to_bed(coordinates, lengths, "g4.motifs.bed")
