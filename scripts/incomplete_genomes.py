"""
incomplete_genomes.py
---------------------
Reads data/combined_metadata.tsv and exports a TSV of all Madagascar isolates
that have missing epidemiological metadata fields (subtype, host, location,
collection date, etc.) OR that are not phylo-ready (missing/short segments).

Output: data/madagascar_missing_metadata.tsv

Run:
    conda activate madagascar_influenza
    python scripts/incomplete_genomes.py
"""

from pathlib import Path
import pandas as pd

DATA_DIR = Path("data")

# Epidemiological fields that should ideally be filled for every isolate
META_FIELDS = {
    "Subtype":          "Subtype (e.g. H5N1)",
    "Host":             "Host species",
    "Country":          "Country",
    "Region":           "Sub-national region",
    "Collection_Date":  "Collection date",
    "Passage_History":  "Passage history (clinical/lab passage)",
    "Lineage":          "Lineage",
    "Clade":            "Clade",
    "Pathogenicity":    "Pathogenicity annotation",
}

df = pd.read_csv(DATA_DIR / "combined_metadata.tsv", sep="\t", low_memory=False)
mdg = df[df["Madagascar"]].copy()

# ── Check each metadata field ─────────────────────────────────────────────────
print(f"\nMadagascar isolates: {len(mdg):,}\n")
print(f"{'Field':<25} {'Missing':>8} {'%':>6}")
print("-" * 42)
for col, label in META_FIELDS.items():
    if col not in mdg.columns:
        continue
    n_missing = mdg[col].isna().sum()
    if col == "Collection_Date":
        # Also flag year-only dates (Month will be NaN) as partial
        n_missing = mdg["Month"].isna().sum()
        label += " (full date)"
    pct = n_missing / len(mdg) * 100
    print(f"{label:<25} {n_missing:>8,} {pct:>6.1f}%")

# ── Build per-isolate missing fields list ────────────────────────────────────
def missing_fields(row):
    missing = []
    for col, label in META_FIELDS.items():
        if col not in row.index:
            continue
        val = row[col]
        if col == "Collection_Date":
            if pd.isna(row.get("Month")):
                missing.append("Collection_Date (partial)")
        elif pd.isna(val) or str(val).strip() == "":
            missing.append(col)
    return ";".join(missing)

mdg["Missing_metadata_fields"] = mdg.apply(missing_fields, axis=1)
has_missing = mdg[mdg["Missing_metadata_fields"] != ""].copy()

print(f"\nIsolates with ≥1 missing metadata field: {len(has_missing):,} / {len(mdg):,}")

# ── Export ────────────────────────────────────────────────────────────────────
out_cols = [
    "Isolate_Id", "Isolate_Name", "Subtype", "Host",
    "Country", "Region", "Collection_Date",
    "Passage_History", "Lineage", "Clade", "Pathogenicity",
    "Segments_in_metadata", "n_segs_sequenced", "phylo_ready",
    "Missing_metadata_fields",
]
out = has_missing[[c for c in out_cols if c in has_missing.columns]].copy()
out = out.sort_values("Missing_metadata_fields")
out_path = DATA_DIR / "madagascar_missing_metadata.tsv"
out.to_csv(out_path, sep="\t", index=False)
print(f"Saved: {out_path}  ({len(out):,} rows)")

