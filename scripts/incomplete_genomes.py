"""
incomplete_genomes.py
---------------------
Reads data/combined_metadata.tsv and exports a TSV of all Madagascar isolates
that are NOT phylo-ready, with columns describing exactly what is missing.

Columns in output:
  Missing_metadata_segs  : segments with no accession ID in GISAID metadata
  Missing_sequence_segs  : segments with no sequence in the FASTA download
  Short_sequence_segs    : segments whose sequence is shorter than the
                           minimum length required for phylogenetics

Run:
    conda activate madagascar_influenza
    python scripts/incomplete_genomes.py
"""

from collections import Counter
from pathlib import Path
import pandas as pd

DATA_DIR = Path("data")
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
MIN_SEQ_LEN = {
    "PB2": 1800, "PB1": 1800, "PA": 1700,
    "HA":  1300, "NP":  1200, "NA": 1000,
    "MP":   800, "NS":   700,
}

df = pd.read_csv(DATA_DIR / "combined_metadata.tsv", sep="\t", low_memory=False)
mdg = df[df["Madagascar"]].copy()
incomplete = mdg[~mdg["phylo_ready"]].copy()


def miss_meta(row):
    return [s for s in SEGMENTS if pd.isna(row.get(f"{s} Segment_Id"))]

def miss_seq(row):
    return [s for s in SEGMENTS if pd.isna(row.get(f"{s}_length"))]

def short_seq(row):
    return [
        s for s in SEGMENTS
        if pd.notna(row.get(f"{s}_length")) and row[f"{s}_length"] < MIN_SEQ_LEN[s]
    ]


incomplete["Missing_metadata_segs"] = incomplete.apply(miss_meta, axis=1)
incomplete["Missing_sequence_segs"] = incomplete.apply(miss_seq,  axis=1)
incomplete["Short_sequence_segs"]   = incomplete.apply(short_seq, axis=1)

# ── Summary ──────────────────────────────────────────────────────────────────
print(f"\nMadagascar isolates : {len(mdg):,}")
print(f"Phylo-ready         : {mdg['phylo_ready'].sum():,}")
print(f"Incomplete          : {len(incomplete):,}")

print("\nMissing from GISAID metadata (segment ID absent):")
for s, c in Counter(s for lst in incomplete["Missing_metadata_segs"] for s in lst).most_common():
    print(f"  {s}: {c:,}")

print("\nAbsent from FASTA (no sequence downloaded):")
for s, c in Counter(s for lst in incomplete["Missing_sequence_segs"] for s in lst).most_common():
    print(f"  {s}: {c:,}")

print("\nSequence present but too short for phylogenetics:")
for s, c in Counter(s for lst in incomplete["Short_sequence_segs"] for s in lst).most_common():
    print(f"  {s}: {c:,}")

# ── Export ────────────────────────────────────────────────────────────────────
out_cols = [
    "Isolate_Id", "Isolate_Name", "Subtype", "Host",
    "Country", "Region", "Collection_Date",
    "Segments_in_metadata", "n_segs_sequenced",
    "Missing_metadata_segs", "Missing_sequence_segs", "Short_sequence_segs",
]
out = incomplete[out_cols].copy()
for col in ["Missing_metadata_segs", "Missing_sequence_segs", "Short_sequence_segs"]:
    out[col] = out[col].apply(lambda x: ";".join(x) if x else "")

out = out.sort_values(["Segments_in_metadata", "n_segs_sequenced"], ascending=False)
out_path = DATA_DIR / "madagascar_incomplete_genomes.tsv"
out.to_csv(out_path, sep="\t", index=False)
print(f"\nSaved: {out_path}  ({len(out):,} rows)")
