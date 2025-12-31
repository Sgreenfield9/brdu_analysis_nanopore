import os
from BCBio import GFF
from dotenv import dotenv_values

def iter_features(features):
    """
    Recursively iterate through all features and their sub-features in a GFF record.

    We will use yield to help save memory since we are looking through a large GFF record.
    Yield will extract each feature one-by-one instead of trying to extracting it all.
    We are also checking for sub_features within the TE's and will recursively
    iterate through all of them 
    """
    for feature in features:
        # We use yield to save memory (used with big .GFF files)
        yield feature
        # If sub-feature exists, recursively iterate through them
        sub_features = getattr(feature, "sub_features", None)
        if sub_features:
            yield from iter_features(sub_features)

# Constant used for TE detection
TE_FEATURE_TYPES = {"mobile_genetic_element", "long_terminal_repeat"}

def get_first_qual(feature, key, default=""):
    """
    Return the first qualifire value as a string, if there is one.
    We have two different classes(qualifiers)
     - Class 1 (Retrotransposons)
       - This contains the Long Terminal Repeats(LTR) Retrotransposons
         and non LTR Retrotransposons
     - Class 2 (DNA Transposons)
    """
    q = feature.qualifiers or {}
    val = q.get(key)
    if not val:
        return default
    if isinstance(val, (list, tuple)):
        return str(val[0]) if val else default
    return str(val)

def is_te_feature(feature, include_ltrs=True):
    """
    Identify TE-related features
    If it is true we will include both TE bodies + LTR features
    If it is false we will include only TE bodies
    """
    ftype = (feature.type or "").strip()

    # Used to find LTR features for the primary rule 
    if include_ltrs:
        if ftype in TE_FEATURE_TYPES:
            return True
    else:
        if ftype == "mobile_genetic_element":
            return True 

    # Secondary rule based on qualifiers present in the file      
    gbkey = get_first_qual(feature, "gbkey", "").strip()

    if gbkey == "mobile_element":
        return True # TE Body
    if include_ltrs and gbkey == "repeat region":
        rpt_type = get_first_qual(feature, "rpt_type", "").strip()
        if rpt_type == "long_terminal_repeat":
            return True
        
    # Signal for TE bodies
    if get_first_qual(feature, "mobile_element_type", ""):
        return True
    
    return False

def get_te_name(feature):
    """
    Naming the TE mobile_element_type -> Name -> ID
    """
    for key in ("mobile_element_type", "Name", "ID"):
        val = get_first_qual(feature, key, "")
        if val:
            return val
    return "TE"

def to_bed_strand(strand_val):
    """Convert Biopython strand to BED strand."""
    return "+" if strand_val == 1 else "-" if strand_val == -1 else "."


def parse_te_gff_to_bed(gff_path, bed_path, include_ltrs=True):
    """
    Parse GFF3 and write TE coordinates to BED.

    Output columns:
      1 chrom
      2 start (0-based)
      3 end (exclusive)
      4 name
      5 length
      6 strand
      7 feature_type
    """
    count = 0
    with open(gff_path) as in_handle, open(bed_path, "w") as out_handle:
        for record in GFF.parse(in_handle):
            chrom = record.id
            for feature in iter_features(record.features):
                if not is_te_feature(feature, include_ltrs=include_ltrs):
                    continue

                start = int(feature.location.start)
                end = int(feature.location.end)
                length = end - start

                name = get_te_name(feature)
                strand = to_bed_strand(feature.location.strand)
                ftype = feature.type or "."

                out_handle.write(
                    f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\t{ftype}\n"
                )
                count += 1

    return count


if __name__ == "__main__":
    # --- Load paths from env/.env (same pattern you used for tRNA) ---
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    env_path = os.path.join(repo_root, "env", ".env")
    env_vars = dotenv_values(env_path)

    gff_file = env_vars.get("W303_GFF")
    if not gff_file:
        raise ValueError("W303_GFF not found in env/.env")

    # Choose output(s)
    # te_bed_all = "w303_te_and_ltrs.bed"
    # te_bed_bodies_only = "w303_te_bodies_only.bed"

    # (A) TE bodies + LTRs (recommended if you want everything TE-related)
    # n_all = parse_te_gff_to_bed(gff_file, te_bed_all, include_ltrs=True)
    # print(f"BED file created: {te_bed_all} ({n_all} features)")

    # (B) TE bodies only (recommended if you only want whole elements, not LTR subparts)
    # n_bodies = parse_te_gff_to_bed(gff_file, te_bed_bodies_only, include_ltrs=False)
    # print(f"BED file created: {te_bed_bodies_only} ({n_bodies} features)")


    