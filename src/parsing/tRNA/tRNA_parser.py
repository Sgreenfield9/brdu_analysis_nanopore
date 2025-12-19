import os

from BCBio import GFF
from dotenv import dotenv_values


def iter_features(features):
    for feature in features:
        yield feature
        sub_features = getattr(feature, "sub_features", None)
        if sub_features:
            for sub_feature in iter_features(sub_features):
                yield sub_feature


repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
env_path = os.path.join(repo_root, "env", ".env")
env_vars = dotenv_values(env_path)
gff_file = env_vars.get("W303_GFF")
if not gff_file:
    raise ValueError("W303_GFF not found in env/.env")

bed_file = "trna_coordinates.bed"

with open(gff_file) as in_handle, open(bed_file, "w") as out_handle:
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
