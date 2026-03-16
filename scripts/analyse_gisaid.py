"""analyse_gisaid.py
------------------
Loads all GISAID EpiFlu metadata exports and FASTA sequence files from
the data/ directory, merges them into a single clean TSV, and prints
a summary of genome completeness.

File discovery rules
--------------------
Metadata  : data/gisaid_epiflu_isolates_<tag>.xls   (one or more)
Sequences : data/gisaid_epiflu_sequence_<tag>.fasta  (zero or more)

The two file types are discovered INDEPENDENTLY, so FASTA files do not
need to have a matching XLS basename. This allows you to split a large
download across multiple FASTA files (e.g. a_excepth3_human_africa.fasta
+ a_h3_human_africa.fasta both contributing sequences to a_human_africa.xls).

Pipeline steps
--------------
1. Load every XLS → parse dates, split Location, compute metadata completeness
2. Keep only Africa records and remove intra-file Isolate_Id duplicates
3. Concatenate all XLS frames → remove cross-file Isolate_Id duplicates
4. Parse every FASTA with seqkit fx2tab → build EPI_ID → length lookup
5. Join sequence lengths onto metadata using the per-segment EPI accessions
   stored in the "<SEG> Segment_Id" columns
6. Compute sequence-level completeness flags
7. Write data/combined_metadata.tsv

Key output columns
------------------
Complete_metadata  : True if all 8 segment accession IDs are present in GISAID
                     (the submitter declared 8 segments, but sequences may be
                     partial or absent from the FASTA download)
<SEG>_length       : nucleotide length of the segment sequence from the FASTA;
                     NaN if no sequence was downloaded for that segment
n_segs_sequenced   : number of core segments (0-8) that have any sequence
complete_sequence  : True if all 8 core segments have a sequence in the FASTA
phylo_ready        : True if all 8 core segments have a sequence AND each
                     exceeds the minimum length threshold (see MIN_SEQ_LEN)

Run
---
    conda activate madagascar_influenza
    python scripts/analyse_gisaid.py
"""

import subprocess
import shutil
import sys
from pathlib import Path

import pandas as pd

# ── Configuration ──────────────────────────────────────────────────────────────

DATA_DIR = Path("data")

# seqkit binary — prefer the conda-env version, fall back to ~/miniconda3
SEQKIT = shutil.which("seqkit") or str(Path.home() / "miniconda3/bin/seqkit")

# The 8 core influenza A segments in polymerase→surface→structural order
CORE_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

# Minimum sequence length (nt) for a segment to be considered usable in
# phylogenetic analysis.  These are set at ~75-80 % of the typical full-length
# reference sequences for influenza A and are used only to flag "phylo_ready".
# They are NOT derived from a single publication; adjust if your analysis
# requires stricter criteria (see README for details).
MIN_SEQ_LEN = {
    "PB2": 1800,   # full length ~2 341 nt
    "PB1": 1800,   # full length ~2 341 nt
    "PA":  1700,   # full length ~2 233 nt
    "HA":  1300,   # full length ~1 700-1 800 nt (varies by subtype)
    "NP":  1200,   # full length ~1 565 nt
    "NA":  1000,   # full length ~1 350-1 470 nt (varies by subtype)
    "MP":   800,   # full length ~1 027 nt
    "NS":   700,   # full length ~  890 nt
}


# ── Metadata loading ──────────────────────────────────────────────────────────

def parse_date(val):
    """Parse the three date formats exported by GISAID EpiFlu.

    GISAID exports dates in one of three formats depending on the precision
    recorded by the submitter:
      - Full date : 'Nov-25-2024'  → format %b-%d-%Y
      - ISO date  : '2024-11-25'   → format %Y-%m-%d
      - Year only : '2024'         → format %Y

    Returns a pandas Timestamp or NaT if parsing fails.
    """
    if pd.isna(val):
        return pd.NaT
    s = str(val).strip()
    for fmt in ("%b-%d-%Y", "%Y-%m-%d", "%Y"):
        try:
            return pd.to_datetime(s, format=fmt)
        except ValueError:
            continue
    # Last resort: let pandas guess (returns NaT if it cannot)
    return pd.to_datetime(s, errors="coerce")


def load_metadata(path: Path) -> pd.DataFrame:
    """Load one GISAID EpiFlu metadata XLS file and return a clean DataFrame.

    Steps performed on each file:
      - Split the 'Location' column ("Continent / Country / Region") into
        separate Continent, Country and Region columns
      - Parse Collection_Date into a proper datetime; extract Year and Month
      - Count how many of the 8 core segment accession IDs are present in
        GISAID metadata (Complete_metadata flag)
      - Drop non-Africa records (all downloads should be Africa-only, but
        this acts as a safety check)
      - Drop intra-file duplicate Isolate_Id rows, keeping the first occurrence
    """
    tag = path.stem.replace("gisaid_epiflu_isolates_", "")
    print(f"  Loading metadata [{tag}]: {path.name}")
    df = pd.read_excel(path)
    print(f"    → {len(df):,} isolates")

    # ── Parse geographic location ──────────────────────────────────────────
    # GISAID stores location as a slash-separated string, e.g.
    # "Africa / Madagascar / Antananarivo"
    loc = df["Location"].str.split(" / ", expand=True)
    df["Continent"] = loc[0].str.strip()
    df["Country"]   = loc[1].str.strip() if 1 in loc.columns else pd.NA
    df["Region"]    = loc[2].str.strip() if 2 in loc.columns else pd.NA

    # ── Parse collection dates ─────────────────────────────────────────────
    df["Collection_Date"] = df["Collection_Date"].apply(parse_date)
    df["Year"]  = df["Collection_Date"].dt.year
    df["Month"] = df["Collection_Date"].dt.month

    # ── Metadata completeness ──────────────────────────────────────────────
    # For each of the 8 core segments, GISAID exports a "<SEG> Segment_Id"
    # column containing the EPI accession of that segment (e.g.
    # "EPI647532|A/goose/Egypt/...").  A non-null value means the submitter
    # deposited at least the segment metadata.  This does NOT guarantee that
    # the actual sequence is in the FASTA download.
    seg_id_cols = [
        f"{s} Segment_Id" for s in CORE_SEGMENTS
        if f"{s} Segment_Id" in df.columns
    ]
    df["Segments_in_metadata"] = df[seg_id_cols].notna().sum(axis=1)
    df["Complete_metadata"]    = df["Segments_in_metadata"] == 8

    # ── Tidy subtype ───────────────────────────────────────────────────────
    # Replace the string 'nan' (which appears when the cell is empty after
    # Excel import) with a proper missing value
    df["Subtype"] = df["Subtype"].astype(str).str.strip().replace("nan", pd.NA)

    df["Source"]     = tag          # which XLS file this row came from
    df["Madagascar"] = df["Country"] == "Madagascar"

    # ── Filter to Africa only ──────────────────────────────────────────────
    # All downloads are filtered to Africa in GISAID, but we double-check
    # here to catch any unexpected records.
    before = len(df)
    df = df[df["Continent"] == "Africa"].copy()
    removed = before - len(df)
    if removed:
        print(f"    ✗ Dropped {removed:,} non-Africa isolates (kept {len(df):,})")
    else:
        print(f"    ✓ All isolates are from Africa")

    # ── Remove intra-file duplicates ───────────────────────────────────────
    # Isolate_Id is GISAID's unique identifier for an influenza isolate.
    # Duplicates within a single export are rare but can happen.
    dups = df.duplicated(subset="Isolate_Id", keep=False)
    if dups.any():
        n_dup        = dups.sum()
        n_unique_dup = df.loc[dups, "Isolate_Id"].nunique()
        print(f"    ⚠ {n_dup:,} rows share {n_unique_dup:,} duplicate Isolate_Id(s) — keeping first occurrence")
        df = df.drop_duplicates(subset="Isolate_Id", keep="first")
        print(f"    → {len(df):,} unique isolates after deduplication")
    else:
        print(f"    ✓ No duplicate Isolate_Id within file")

    return df


# ── FASTA parsing ─────────────────────────────────────────────────────────────

def parse_fasta_lengths(fasta_path: Path) -> pd.DataFrame:
    """Extract EPI accession, segment name, and sequence length from a FASTA.

    Uses seqkit fx2tab (fast C binary) instead of pure-Python parsing because
    the FASTA files can be very large (>100 MB).

    GISAID EpiFlu FASTA headers follow this format:
      >EPI3754149|NA|A/chicken/Mozambique/857-P3-2/2023|EPI_ISL_19644683|A_/_H7N6
       ^^^^^^^^^  ^^
       EPI_ID     Segment name (field [1])

    Returns a DataFrame with columns: EPI_ID, Segment, Length.
    """
    print(f"  Parsing FASTA [{fasta_path.name}] with seqkit ...")
    result = subprocess.run(
        [SEQKIT, "fx2tab", "--name", "--length", str(fasta_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"    ✗ seqkit error: {result.stderr[:200]}", file=sys.stderr)
        return pd.DataFrame(columns=["EPI_ID", "Segment", "Length"])

    rows = []
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        # seqkit fx2tab outputs: <header>\t<length>
        header, length = line.rsplit("\t", 1)
        fields = header.split("|")
        if len(fields) >= 2:
            rows.append({
                "EPI_ID":  fields[0].strip(),   # bare EPI accession
                "Segment": fields[1].strip(),   # segment name (PB2, HA, etc.)
                "Length":  int(length),
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        print(f"    → {len(df):,} sequences, {df['Segment'].nunique()} segment types")
    return df


def attach_sequence_lengths(meta: pd.DataFrame, fasta_df: pd.DataFrame) -> pd.DataFrame:
    """Join per-segment sequence lengths onto the metadata DataFrame.

    How the join works
    ------------------
    Each core segment has a dedicated column in the GISAID metadata:
      "HA Segment_Id" → 'EPI647532|A/goose/Egypt/BSU-NLQP-DAK-2/2009'

    The EPI accession before the pipe ('EPI647532') is the unique GISAID
    identifier for that specific sequence.  The FASTA headers use the same
    bare EPI accession as their first pipe-delimited field:
      >EPI647532|HA|A/goose/Egypt/...

    We therefore:
      1. Extract the bare accession from each '<SEG> Segment_Id' column
      2. Look it up in the fasta_df EPI_ID → Length mapping
      3. Store the result in '<SEG>_length'

    Isolates without a matching FASTA sequence get NaN for that segment's length.
    """
    if fasta_df.empty:
        # No sequences available at all; initialise all length columns to NaN
        for seg in CORE_SEGMENTS:
            meta[f"{seg}_length"] = pd.NA
        return meta

    # Build a Series: EPI_ID (index) → Length (value)
    # Each EPI accession is globally unique in GISAID, so no duplicate-key issue.
    epi_to_len = fasta_df.set_index("EPI_ID")["Length"]

    for seg in CORE_SEGMENTS:
        col_id  = f"{seg} Segment_Id"   # metadata column with EPI accession
        col_len = f"{seg}_length"        # new column we are adding
        if col_id in meta.columns:
            # Strip everything after the pipe to get the bare EPI accession
            bare = meta[col_id].str.split("|").str[0]
            meta[col_len] = bare.map(epi_to_len)
        else:
            # Segment ID column not present in this export → all NaN
            meta[col_len] = pd.NA

    return meta


def add_sequence_completeness(df: pd.DataFrame) -> pd.DataFrame:
    """Compute sequence-level completeness columns.

    Adds three columns:
      n_segs_sequenced  : int 0-8, how many core segments have any sequence
      complete_sequence : True if all 8 core segments have a sequence
      phylo_ready       : True if all 8 core segments have a sequence AND
                          each sequence meets the MIN_SEQ_LEN threshold

    Note: complete_sequence only checks presence (length > 0), while
    phylo_ready additionally enforces minimum lengths suitable for
    building maximum-likelihood phylogenetic trees.
    """
    len_cols = [f"{s}_length" for s in CORE_SEGMENTS if f"{s}_length" in df.columns]

    # Count how many segments have any sequence (non-NaN length)
    df["n_segs_sequenced"]  = df[len_cols].notna().sum(axis=1)
    df["complete_sequence"] = df["n_segs_sequenced"] == 8

    # phylo_ready: every segment present AND above its minimum length
    df["phylo_ready"] = df[len_cols].apply(
        lambda row: all(
            pd.notna(row[f"{seg}_length"]) and row[f"{seg}_length"] >= MIN_SEQ_LEN[seg]
            for seg in CORE_SEGMENTS
            if f"{seg}_length" in df.columns
        ),
        axis=1,
    )
    return df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # ── Step 1: Discover and load all metadata XLS files ──────────────────
    xls_files = sorted(DATA_DIR.glob("gisaid_epiflu_isolates_*.xls"))
    if not xls_files:
        sys.exit("No gisaid_epiflu_isolates_*.xls files found in data/")

    print(f"\nFound {len(xls_files)} metadata file(s):")
    frames       = [load_metadata(p) for p in xls_files]
    total_before = sum(len(f) for f in frames)
    combined     = pd.concat(frames, ignore_index=True)

    # ── Step 2: Remove cross-file duplicates ──────────────────────────────
    # The same isolate can appear in multiple GISAID downloads (e.g. it is
    # both a human isolate AND features in the animal download).  We keep
    # the first occurrence, which corresponds to the alphabetically first XLS.
    cross_dups = combined.duplicated(subset="Isolate_Id", keep=False)
    if cross_dups.any():
        n_dup        = cross_dups.sum()
        n_unique_dup = combined.loc[cross_dups, "Isolate_Id"].nunique()
        print(f"\n⚠ Cross-file duplicates: {n_dup:,} rows share {n_unique_dup:,} Isolate_Id(s) — keeping first occurrence")
        combined = combined.drop_duplicates(subset="Isolate_Id", keep="first")
    else:
        print(f"\n✓ No cross-file duplicate Isolate_Id(s)")
    print(f"Combined (deduplicated): {len(combined):,} isolates (from {total_before:,} across all files)")

    # ── Step 3: Discover and parse ALL FASTA files ────────────────────────
    # FASTA files are discovered INDEPENDENTLY of the XLS files.  This means
    # a single XLS (e.g. a_human_africa.xls) can have its sequences split
    # across multiple FASTAs (e.g. a_excepth3_human_africa.fasta +
    # a_h3_human_africa.fasta) without any naming-convention requirements.
    fasta_files = sorted(DATA_DIR.glob("gisaid_epiflu_sequence_*.fasta"))
    xls_tags    = {p.stem.replace("gisaid_epiflu_isolates_", "") for p in xls_files}

    all_seqs = []
    print()
    for fasta_path in fasta_files:
        tag = fasta_path.stem.replace("gisaid_epiflu_sequence_", "")
        if tag not in xls_tags:
            # FASTA has no direct matching XLS — still parse it; its EPI
            # accessions will match isolates already loaded from other XLS files
            print(f"  Note: [{fasta_path.name}] has no matching XLS — sequences will be joined by EPI accession")
        all_seqs.append(parse_fasta_lengths(fasta_path))

    if all_seqs:
        # Merge all FASTA records; if the same EPI_ID appears in multiple
        # FASTAs (batch split), keep the first occurrence
        fasta_df = pd.concat(all_seqs, ignore_index=True).drop_duplicates("EPI_ID")
        print(f"\n  Total unique sequences across all FASTAs: {len(fasta_df):,}")
        combined = attach_sequence_lengths(combined, fasta_df)
    else:
        print("  No FASTA files found — sequence length columns will be NaN")
        for seg in CORE_SEGMENTS:
            combined[f"{seg}_length"] = pd.NA

    # ── Step 4: Compute completeness flags ────────────────────────────────
    combined = add_sequence_completeness(combined)

    # ── Summary ───────────────────────────────────────────────────────────────
    africa = combined[combined["Continent"] == "Africa"]
    mdg    = combined[combined["Madagascar"]]

    for label, sub in [("Madagascar", mdg), ("All Africa", africa)]:
        n  = len(sub)
        cm = int(sub["Complete_metadata"].sum())
        cs = int(sub["complete_sequence"].sum()) if "complete_sequence" in sub else 0
        pr = int(sub["phylo_ready"].sum()) if "phylo_ready" in sub else 0
        yr = int(sub["Year"].between(1990, 2026).sum())
        print(f"\n── {label} (n={n:,}) ──")
        print(f"  Complete metadata (8 seg IDs) : {cm:,} ({cm/n*100:.1f}%)")
        print(f"  Complete sequences (8 seqs)   : {cs:,} ({cs/n*100:.1f}%)")
        print(f"  Phylo-ready (min lengths)     : {pr:,} ({pr/n*100:.1f}%)")
        print(f"  Valid collection date          : {yr:,} ({yr/n*100:.1f}%)")

    # ── Save ──────────────────────────────────────────────────────────────────
    out = DATA_DIR / "combined_metadata.tsv"
    combined.to_csv(out, sep="\t", index=False)
    print(f"\nSaved: {out}  ({len(combined):,} rows, {len(combined.columns)} columns)")


if __name__ == "__main__":
    main()
