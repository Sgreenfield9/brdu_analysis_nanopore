import pysam
import argparse
import numpy as np


def phred_to_prob(phred_score):
    """Convert phred score to probability of modification"""
    return 1 - 10**(-phred_score / 10)


def parse_dnascent_mm_ml(mm_tag, ml_tag, mod_type='b'):
    """
    Parse DNAscent MM and ML tags to extract modification data for specific type
    
    Args:
        mm_tag: MM tag string (e.g., "N+b?,1,2,3;N+e?,4,5,6")
        ml_tag: ML tag array
        mod_type: 'b' for BrdU or 'e' for EdU
    
    Returns:
        (skips, scores) tuple or (None, None) if not found
    """
    # Split by semicolon to handle multiple modification types
    mod_sections = mm_tag.split(';')
    
    target_prefix = f"N+{mod_type}?"
    
    # Track which section we want and where it is
    section_idx = -1
    skips = None
    
    for idx, section in enumerate(mod_sections):
        if section.strip().startswith(target_prefix):
            section_idx = idx
            parts = section.split(',')
            if len(parts) < 2:
                return None, None
            try:
                skips = [int(x) for x in parts[1:]]
            except ValueError:
                return None, None
            break
    
    if skips is None:
        return None, None
    
    # ML tag contains scores for all modification types concatenated
    # We need to figure out which scores belong to our modification type
    # Parse all sections to know how many scores each has
    all_skip_counts = []
    for section in mod_sections:
        if not section.strip():
            continue
        parts = section.split(',')
        if len(parts) >= 2:
            try:
                section_skips = [int(x) for x in parts[1:]]
                all_skip_counts.append(len(section_skips))
            except ValueError:
                all_skip_counts.append(0)
    
    # Calculate offset into ML array
    ml_offset = sum(all_skip_counts[:section_idx])
    ml_length = len(skips)
    
    # Extract the relevant ML scores
    if ml_offset + ml_length <= len(ml_tag):
        scores = ml_tag[ml_offset:ml_offset + ml_length]
    else:
        scores = ml_tag[ml_offset:]
    
    return skips, scores


def count_mod(bam_file, probability_threshold, mod_type='b'):
    """
    Count number of T bases with modification probability >= threshold
    
    Args:
        bam_file: Path to mod.bam file
        probability_threshold: Minimum probability threshold (0-1)
        mod_type: 'b' for BrdU or 'e' for EdU
    
    Returns:
        Dictionary with read_id as key and count as value
    """
    # Convert probability to phred score threshold
    if probability_threshold >= 1.0:
        phred_threshold = 255
    elif probability_threshold <= 0:
        phred_threshold = 0
    else:
        phred_threshold = -10 * np.log10(1 - probability_threshold)
    
    results = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            read_id = read.query_name
            
            # Get modification tags
            if not read.has_tag("MM") or not read.has_tag("ML"):
                results[read_id] = 0
                continue
            
            mm_tag = read.get_tag("MM")
            ml_tag = read.get_tag("ML")
            
            # Parse DNAscent format
            skips, scores = parse_dnascent_mm_ml(mm_tag, ml_tag, mod_type)
            
            if skips is None or scores is None:
                results[read_id] = 0
                continue
            
            # Count modifications above threshold
            count = sum(1 for score in scores if score >= phred_threshold)
            results[read_id] = count
    
    return results


def sliding_window_mod_to_bed6(bam_file, probability_threshold, window_size, mod_type='b'):
    """
    Calculate smoothed modification using sliding window and output as BED6
    
    Args:
        bam_file: Path to mod.bam file
        probability_threshold: Minimum probability threshold (0-1)
        window_size: Number of T bases in sliding window
        mod_type: 'b' for BrdU or 'e' for EdU
    
    Returns:
        List of BED6 tuples: (chrom, start, end, name, score, strand)
    """
    # Convert probability to phred score threshold
    if probability_threshold >= 1.0:
        phred_threshold = 255
    elif probability_threshold <= 0:
        phred_threshold = 0
    else:
        phred_threshold = -10 * np.log10(1 - probability_threshold)
    
    bed_entries = []
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            read_id = read.query_name
            chrom = read.reference_name
            
            if chrom is None:
                continue
            
            # Get strand
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            
            # Get modification tags
            if not read.has_tag("MM") or not read.has_tag("ML"):
                continue
            
            mm_tag = read.get_tag("MM")
            ml_tag = read.get_tag("ML")
            
            # Parse DNAscent format
            skips, scores = parse_dnascent_mm_ml(mm_tag, ml_tag, mod_type)
            
            if skips is None or scores is None:
                continue
            
            # Find T positions in the sequence
            sequence = read.query_sequence
            if sequence is None:
                continue
            
            t_positions = [i for i, base in enumerate(sequence) if base == 'T']
            
            # Map ML scores to T positions using skip counts
            mod_data = []  # (query_pos, ref_pos, score)
            current_t_idx = 0
            
            # Compute reference positions ONCE per read
            ref_positions = read.get_reference_positions()

            for skip, score in zip(skips, scores):
                current_t_idx += skip
                if current_t_idx < len(t_positions):
                    query_pos = t_positions[current_t_idx]
                    ref_pos = ref_positions[query_pos] if query_pos < len(ref_positions) else None
                    if ref_pos is not None:
                        mod_data.append((query_pos, ref_pos, score))
                    current_t_idx += 1
            
            if len(mod_data) == 0:
                continue
            
            # Apply sliding window
            for i in range(len(mod_data)):
                # Get window of scores
                window_start = max(0, i - window_size // 2)
                window_end = min(len(mod_data), i + window_size // 2 + 1)
                window_scores = [mod_data[j][2] for j in range(window_start, window_end)]
                
                # Calculate proportion above threshold in window
                above_threshold = sum(1 for score in window_scores if score >= phred_threshold)
                smoothed_prob = above_threshold / len(window_scores)
                
                # Get reference position and window boundaries
                ref_pos = mod_data[i][1]
                window_ref_start = mod_data[window_start][1]
                window_ref_end = mod_data[window_end - 1][1] + 1  # BED is half-open
                
                # BED6 format: chrom, start, end, name, score, strand
                # Score should be 0-1000 for BED format
                bed_score = int(smoothed_prob * 1000)
                bed_name = f"{read_id}_T{i}"
                
                bed_entries.append((chrom, window_ref_start, window_ref_end, 
                                   bed_name, bed_score, strand))
    
    return bed_entries


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze BrdU/EdU modifications from DNAscent mod.bam file")
    parser.add_argument("--bam", required=True, help="Path to DNAscent mod.bam file")
    parser.add_argument("--prob", type=float, required=True, 
                        help="Probability threshold (0-1). DNAscent recommends 0.5")
    parser.add_argument("--mod", choices=['b', 'e', 'brdu', 'edu'], default='b',
                        help="Modification type: 'b'/'brdu' for BrdU, 'e'/'edu' for EdU (default: b)")
    parser.add_argument("--window", type=int, default=None,
                        help="Window size for sliding window analysis (number of T bases)")
    parser.add_argument("--output", default=None,
                        help="Output file (default: print to stdout)")
    parser.add_argument("--format", choices=['bed6', 'tsv'], default='tsv',
                        help="Output format: bed6 or tsv (default: tsv)")
    
    args = parser.parse_args()
    
    # Validate probability
    if not 0 <= args.prob <= 1:
        raise ValueError("Probability must be between 0 and 1")
    
    # Normalize mod type
    mod_type = 'b' if args.mod in ['b', 'brdu'] else 'e'
    mod_name = 'BrdU' if mod_type == 'b' else 'EdU'
    
    if args.window is None:
        # Simple count mode
        print(f"Counting T bases with {mod_name} probability >= {args.prob}")
        results = count_mod(args.bam, args.prob, mod_type)
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write("read_id\tcount\n")
                for read_id, count in results.items():
                    f.write(f"{read_id}\t{count}\n")
        else:
            print("read_id\tcount")
            for read_id, count in results.items():
                print(f"{read_id}\t{count}")
    
    else:
        # Sliding window mode
        print(f"Calculating smoothed {mod_name} with window size {args.window} and threshold {args.prob}")
        bed_entries = sliding_window_mod_to_bed6(args.bam, args.prob, args.window, mod_type)
        
        if args.format == 'bed6':
            # Output as BED6
            if args.output:
                with open(args.output, 'w') as f:
                    for entry in bed_entries:
                        f.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}\n")
            else:
                for entry in bed_entries:
                    print(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}")
        else:
            # Output as TSV with header
            if args.output:
                with open(args.output, 'w') as f:
                    f.write("chrom\tstart\tend\tname\tscore\tstrand\n")
                    for entry in bed_entries:
                        f.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}\n")
            else:
                print("chrom\tstart\tend\tname\tscore\tstrand")
                for entry in bed_entries:
                    print(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}")