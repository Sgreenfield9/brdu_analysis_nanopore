import pandas as pd
import matplotlib.pyplot as plt

# ==========================
# User inputs
# ==========================
bed_file = "/Users/samuelgreenfiled/Desktop/SG_EXP_2025_010/sort_index/S_phase_brdu_smoothed.bed"

# Optional: zoom into a region (chromosomal coordinates)
# Set to None to plot full contigs
zoom_region = None  # e.g., (14500, 15000)
# ==========================

# Load BED6 file
bed = pd.read_csv(bed_file, sep='\t', header=None,
                  names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

# Convert score to probability (0-1)
bed['prob'] = bed['score'] / 1000.0

# Loop over each contig
for chrom in bed['chrom'].unique():
    chrom_bed = bed[bed['chrom'] == chrom].copy()

    # Apply zoom if specified
    if zoom_region:
        start, end = zoom_region
        chrom_bed = chrom_bed[(chrom_bed['start'] >= start) & (chrom_bed['end'] <= end)]
        if chrom_bed.empty:
            continue  # skip if nothing in region

    # Split by strand
    plus_strand = chrom_bed[chrom_bed['strand'] == '+']
    minus_strand = chrom_bed[chrom_bed['strand'] == '-']

    # Create figure
    plt.figure(figsize=(15, 3))

    # Plot + strand
    if not plus_strand.empty:
        plt.plot(plus_strand['start'], plus_strand['prob'], color='green', linewidth=1)
        plt.fill_between(plus_strand['start'], 0, plus_strand['prob'], color='green', alpha=0.3)

    # Plot - strand
    if not minus_strand.empty:
        plt.plot(minus_strand['start'], minus_strand['prob'], color='red', linewidth=1)
        plt.fill_between(minus_strand['start'], 0, minus_strand['prob'], color='red', alpha=0.3)

    plt.xlabel('Genomic position')
    plt.ylabel('BrdU probability (sliding window)')
    plt.ylim(0, 1)
    plt.title(f'Smoothed BrdU signal along {chrom} (+ green / - red)')

    plt.tight_layout()
    plt.show()
