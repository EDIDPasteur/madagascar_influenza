"""analyse_gisaid.py
------------------
Loads all GISAID EpiFlu metadata exports and FASTA sequence files from
the data/ directory, merges them into a single clean TSV, and writes a
provenance sidecar JSON.

File discovery rules
--------------------
Metadata  : data/gisaid_epiflu_isolates_<tag>.xls   (one or more)
Sequences : data/gisaid_epiflu_sequence_<tag>.fasta  (zero or more, independent
            of XLS names — a single XLS can have its FASTA split across multiple
            files, e.g. a_excepth3_human_africa.fasta + a_h3_human_africa.fasta)

Pipeline steps
--------------
1.  Validate that seqkit is available.
2.  Load every XLS → parse dates, split Location, compute metadata completeness.
3.  Keep only Africa records; remove intra-file Isolate_Id duplicates.
4.  Concatenate all XLS frames; remove cross-file duplicates (keep first,
    where "first" = alphabetical filename order, then original row order).
5.  Parse every FASTA with seqkit fx2tab → EPI_ID → length lookup.
6.  Join sequence lengths onto metadata via <SEG> Segment_Id columns.
7.  Compute sequence-level completeness flags (vectorised).
8.  Sort by Isolate_Id for deterministic output.
9.  Write data/combined_metadata.tsv and data/combined_metadata.provenance.json.

Key output columns
------------------
Complete_metadata  : all 8 segment accession IDs present in GISAID metadata
<SEG>_length       : nucleotide length from FASTA (np.nan if not downloaded)
n_segs_sequenced   : number of core segments (0-8) with any sequence
complete_sequence  : all 8 core segments have a sequence
phylo_ready        : all 8 core segments have a sequence AND each meets MIN_SEQ_LEN

Run
---
    conda activate madagascar_influenza
    python scripts/analyse_gisaid.py [--data-dir data] [--output data/combined_metadata.tsv]
"""

import argparse
import json
import logging
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# -- Logging ------------------------------------------------------------------
# Emit INFO to stderr AND a log file next to the output TSV.
def _setup_logging(log_path: Path) -> None:
    fmt = "%(asctime)s  %(levelname)-8s  %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
        ],
    )

log = logging.getLogger(__name__)

# -- Configuration ------------------------------------------------------------

CORE_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

# Minimum sequence length (nt) for a segment to be considered usable in
# phylogenetic analysis (~75-80% of full-length reference).
# Adjust to match the target pipeline (e.g. Nextstrain avian-flu augur filter).
MIN_SEQ_LEN: dict = {
    "PB2": 1800,  # full ~2341 nt
    "PB1": 1800,  # full ~2341 nt
    "PA":  1700,  # full ~2233 nt
    "HA":  1300,  # full ~1700-1800 nt (varies by subtype)
    "NP":  1200,  # full ~1565 nt
    "NA":  1000,  # full ~1350-1470 nt (varies by subtype)
    "MP":   800,  # full ~1027 nt
    "NS":   700,  # full ~890 nt
}

SEQKIT_TIMEOUT = 600  # seconds; increase for very large FASTA files


# -- seqkit validation --------------------------------------------------------

def resolve_seqkit() -> str:
    """Return path to seqkit binary; abort early with a clear message if absent."""
    candidate = shutil.which("seqkit") or str(Path.home() / "miniconda3/bin/seqkit")
    if not Path(candidate).is_file():
        log.error(
            "seqkit binary not found.\n"
            "  Searched PATH via shutil.which() and %s\n"
            "  Install with:  conda install -c bioconda seqkit",
            candidate,
        )
        sys.exit(1)
    log.info("seqkit resolved -> %s", candidate)
    return candidate


# -- Metadata loading ---------------------------------------------------------

def parse_date(val):
    """Parse the three date formats exported by GISAID EpiFlu.

    Formats tried in order:
      'Nov-25-2024'  -> %b-%d-%Y  (full date)
      '2024-11-25'   -> %Y-%m-%d  (ISO date)
      '2024'         -> %Y        (year only; Month will remain NaN)
    Returns pd.Timestamp or pd.NaT.
    """
    if pd.isna(val):
        return pd.NaT
    s = str(val).strip()
    for fmt in ("%b-%d-%Y", "%Y-%m-%d", "%Y"):
        try:
            return pd.to_datetime(s, format=fmt)
        except ValueError:
            continue
    return pd.to_datetime(s, errors="coerce")


def load_metadata(path: Path) -> pd.DataFrame:
    """Load one GISAID EpiFlu metadata XLS and return a clean DataFrame.

    Steps performed:
    - Parse 'Location' string ("Continent / Country / Region") into columns
    - Parse 'Collection_Date' into datetime; extract Year and Month
    - Count segment accession IDs -> Segments_in_metadata, Complete_metadata
    - Filter to Africa only (safety check; downloads should already be Africa)
    - Remove intra-file Isolate_Id duplicates (keep first occurrence)
    """
    tag = path.stem.replace("gisaid_epiflu_isolates_", "")
    log.info("Loading metadata [%s]: %s", tag, path.name)

    try:
        df = pd.read_excel(path)
    except Exception as exc:
        log.error("Failed to read %s: %s -- skipping file.", path.name, exc)
        return pd.DataFrame()

    log.info("  -> %d isolates read", len(df))

    # Geographic location
    # Format: "Africa / Madagascar / Antananarivo"
    if "Location" in df.columns:
        loc = df["Location"].str.split(" / ", expand=True)
        df["Continent"] = loc[0].str.strip()
        df["Country"]   = loc[1].str.strip() if 1 in loc.columns else np.nan
        df["Region"]    = loc[2].str.strip() if 2 in loc.columns else np.nan
    else:
        log.warning("  'Location' column absent -- Continent/Country/Region set to NaN")
        df["Continent"] = df["Country"] = df["Region"] = np.nan

    # Collection date
    if "Collection_Date" in df.columns:
        df["Collection_Date"] = df["Collection_Date"].apply(parse_date)
    else:
        log.warning("  'Collection_Date' column absent -- setting to NaT")
        df["Collection_Date"] = pd.NaT
    df["Year"]  = df["Collection_Date"].dt.year
    df["Month"] = df["Collection_Date"].dt.month

    # Metadata completeness:
    # "<SEG> Segment_Id" is non-null when the submitter deposited that segment.
    # This does NOT guarantee the sequence is present in the FASTA download.
    seg_id_cols = [
        f"{s} Segment_Id" for s in CORE_SEGMENTS
        if f"{s} Segment_Id" in df.columns
    ]
    if not seg_id_cols:
        log.warning(
            "  No '<SEG> Segment_Id' columns found in %s -- "
            "Complete_metadata will always be False. "
            "Check that the XLS was exported with segment columns enabled.",
            path.name,
        )
    df["Segments_in_metadata"] = df[seg_id_cols].notna().sum(axis=1) if seg_id_cols else 0
    df["Complete_metadata"]    = df["Segments_in_metadata"] == 8

    # Subtype: keep genuine NaN rather than converting to string "nan"
    if "Subtype" in df.columns:
        df["Subtype"] = df["Subtype"].where(df["Subtype"].notna(), other=np.nan)
        df["Subtype"] = df["Subtype"].astype(object).str.strip()
    else:
        df["Subtype"] = np.nan

    df["Source"]     = tag
    df["Madagascar"] = df["Country"] == "Madagascar"

    # Africa filter (all downloads should already be Africa-only; this is a
    # safety check to catch any unexpected records)
    before = len(df)
    df = df[df["Continent"] == "Africa"].copy()
    removed = before - len(df)
    if removed:
        log.warning("  Dropped %d non-Africa isolates (kept %d)", removed, len(df))
    else:
        log.info("  All isolates are from Africa")

    # Intra-file deduplication on Isolate_Id (GISAID's unique isolate key)
    dups = df.duplicated(subset="Isolate_Id", keep=False)
    if dups.any():
        n_dup        = int(dups.sum())
        n_unique_dup = int(df.loc[dups, "Isolate_Id"].nunique())
        log.warning(
            "  %d rows share %d Isolate_Id(s) -- keeping first occurrence",
            n_dup, n_unique_dup,
        )
        df = df.drop_duplicates(subset="Isolate_Id", keep="first")
        log.info("  -> %d unique isolates after deduplication", len(df))
    else:
        log.info("  No duplicate Isolate_Id within file")

    return df


# -- FASTA parsing ------------------------------------------------------------

def parse_fasta_lengths(fasta_path: Path, seqkit: str) -> pd.DataFrame:
    """Extract EPI accession, segment name, and length from a GISAID FASTA.

    GISAID EpiFlu FASTA header format:
      >EPI3754149|NA|A/chicken/Mozambique/857-P3-2/2023|EPI_ISL_19644683|A_/_H7N6
       ^^^^^^^^^  ^^
       EPI_ID     Segment name

    Uses seqkit fx2tab for speed on large files (avoids loading sequences into
    Python memory). Returns DataFrame with columns: EPI_ID, Segment, Length.
    """
    log.info("Parsing FASTA [%s] with seqkit ...", fasta_path.name)
    try:
        result = subprocess.run(
            [seqkit, "fx2tab", "--name", "--length", str(fasta_path)],
            capture_output=True,
            text=True,
            timeout=SEQKIT_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        log.error(
            "seqkit timed out after %ds on %s -- skipping FASTA.",
            SEQKIT_TIMEOUT, fasta_path.name,
        )
        return pd.DataFrame(columns=["EPI_ID", "Segment", "Length"])
    except FileNotFoundError:
        log.error("seqkit binary not found at '%s'.", seqkit)
        return pd.DataFrame(columns=["EPI_ID", "Segment", "Length"])

    if result.returncode != 0:
        log.error("seqkit error on %s:\n%s", fasta_path.name, result.stderr[:400])
        return pd.DataFrame(columns=["EPI_ID", "Segment", "Length"])

    rows = []
    skipped = 0
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        if "\t" not in line:
            log.warning("Malformed seqkit output line (no tab): %r -- skipping", line[:80])
            skipped += 1
            continue
        header, length_str = line.rsplit("\t", 1)
        fields = header.split("|")
        if len(fields) < 2:
            log.warning("FASTA header has no pipe delimiter: %r -- skipping", header[:80])
            skipped += 1
            continue
        try:
            rows.append({
                "EPI_ID":  fields[0].strip(),   # bare EPI accession
                "Segment": fields[1].strip(),   # segment name
                "Length":  int(length_str),
            })
        except ValueError:
            log.warning("Could not parse length %r in line %r -- skipping", length_str, line[:80])
            skipped += 1

    if skipped:
        log.warning("Skipped %d malformed lines in %s", skipped, fasta_path.name)

    df = pd.DataFrame(rows)
    if df.empty:
        log.warning("No sequences parsed from %s", fasta_path.name)
    else:
        log.info("  -> %d sequences, %d segment types", len(df), df["Segment"].nunique())
    return df


# -- Sequence length join -----------------------------------------------------

def attach_sequence_lengths(meta: pd.DataFrame, fasta_df: pd.DataFrame) -> pd.DataFrame:
    """Join per-segment sequence lengths onto metadata.

    Join key: bare EPI accession extracted from the "<SEG> Segment_Id" column.
      Metadata value: 'EPI647532|A/goose/Egypt/...'  -> bare key 'EPI647532'
      FASTA EPI_ID:   'EPI647532'

    All length columns are float64 with np.nan for missing values, avoiding
    pd.NA/float mixing issues in downstream numeric comparisons.
    """
    # EPI_ID -> Length (float so NaN propagates correctly)
    epi_to_len = fasta_df.set_index("EPI_ID")["Length"].astype(float)

    for seg in CORE_SEGMENTS:
        col_id  = f"{seg} Segment_Id"
        col_len = f"{seg}_length"
        if col_id in meta.columns:
            # Strip everything after the pipe to get the bare EPI accession
            bare = meta[col_id].str.split("|").str[0]
            meta[col_len] = bare.map(epi_to_len)   # float64, NaN for no match
        else:
            meta[col_len] = np.nan

    return meta


# -- Completeness flags -------------------------------------------------------

def add_sequence_completeness(df: pd.DataFrame) -> pd.DataFrame:
    """Add n_segs_sequenced, complete_sequence, and phylo_ready columns.

    phylo_ready uses vectorised boolean operations (~100x faster than
    row-wise apply and avoids closure-over-df fragility).
    """
    len_cols = [f"{s}_length" for s in CORE_SEGMENTS if f"{s}_length" in df.columns]

    df["n_segs_sequenced"]  = df[len_cols].notna().sum(axis=1).astype(int)
    df["complete_sequence"] = df["n_segs_sequenced"] == len(CORE_SEGMENTS)

    # Vectorised phylo_ready: every segment present AND >= minimum length
    phylo = pd.Series(True, index=df.index)
    for seg in CORE_SEGMENTS:
        col = f"{seg}_length"
        if col in df.columns:
            phylo &= df[col].notna() & (df[col] >= MIN_SEQ_LEN[seg])
        else:
            # Column absent entirely -> can never be phylo_ready
            phylo = pd.Series(False, index=df.index)
            break

    df["phylo_ready"] = phylo
    return df


# -- Provenance ---------------------------------------------------------------

def _seqkit_version(seqkit: str) -> str:
    try:
        r = subprocess.run([seqkit, "version"], capture_output=True, text=True, timeout=10)
        return r.stdout.strip() or r.stderr.strip()
    except Exception:
        return "unknown"


def write_provenance(
    out_path: Path,
    seqkit: str,
    xls_stats: list,
    fasta_stats: list,
    stage_counts: dict,
    min_seq_len: dict,
) -> None:
    """Write a JSON provenance sidecar alongside the output TSV."""
    import platform
    prov = {
        "run_timestamp_utc":     datetime.now(timezone.utc).isoformat(),
        "python_version":        sys.version,
        "pandas_version":        pd.__version__,
        "numpy_version":         np.__version__,
        "seqkit_path":           seqkit,
        "seqkit_version":        _seqkit_version(seqkit),
        "platform":              platform.platform(),
        "min_seq_len_thresholds": min_seq_len,
        "input_xls_files":       xls_stats,
        "input_fasta_files":     fasta_stats,
        "pipeline_stage_counts": stage_counts,
    }
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5)
        prov["git_commit"] = r.stdout.strip() if r.returncode == 0 else "not a git repo"
    except Exception:
        prov["git_commit"] = "unknown"

    prov_path = out_path.with_suffix(".provenance.json")
    with open(prov_path, "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2)
    log.info("Provenance written -> %s", prov_path)


# -- Summary ------------------------------------------------------------------

def _pct(n: int, total: int) -> str:
    if total == 0:
        return "n/a"
    return f"{n:,} ({n / total * 100:.1f}%)"


def print_summary(label: str, sub: pd.DataFrame) -> None:
    n = len(sub)
    if n == 0:
        log.info("-- %s (n=0) -- no records", label)
        return
    cm = int(sub["Complete_metadata"].sum()) if "Complete_metadata" in sub else 0
    cs = int(sub["complete_sequence"].sum()) if "complete_sequence" in sub else 0
    pr = int(sub["phylo_ready"].sum())        if "phylo_ready"       in sub else 0
    yr = int(sub["Year"].between(1900, 2030).sum()) if "Year" in sub else 0
    log.info(
        "\n-- %s (n=%s) --\n"
        "  Complete metadata (8 seg IDs) : %s\n"
        "  Complete sequences (8 seqs)   : %s\n"
        "  Phylo-ready (min lengths)     : %s\n"
        "  Valid collection date         : %s",
        label, f"{n:,}",
        _pct(cm, n), _pct(cs, n), _pct(pr, n), _pct(yr, n),
    )


# -- CLI ----------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merge GISAID EpiFlu metadata XLS + FASTA files into a single TSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--data-dir", type=Path, default=Path("data"),
        help="Directory containing gisaid_epiflu_isolates_*.xls and "
             "gisaid_epiflu_sequence_*.fasta files.",
    )
    p.add_argument(
        "--output", type=Path, default=None,
        help="Output TSV path. Defaults to <data-dir>/combined_metadata.tsv.",
    )
    return p.parse_args()


# -- Main ---------------------------------------------------------------------

def main() -> None:
    args     = parse_args()
    data_dir = args.data_dir
    out_path = args.output or (data_dir / "combined_metadata.tsv")
    log_path = out_path.with_suffix(".log")
    _setup_logging(log_path)

    log.info("=== analyse_gisaid.py ===")
    log.info("Data directory : %s", data_dir.resolve())
    log.info("Output         : %s", out_path.resolve())
    log.info("Log file       : %s", log_path.resolve())

    seqkit = resolve_seqkit()   # exits immediately if not found

    # Step 1: Discover + load all XLS files ----------------------------------
    xls_glob  = "gisaid_epiflu_isolates_*.xls"
    xls_files = sorted(data_dir.glob(xls_glob))
    log.info("XLS glob       : %s/%s  -> %d file(s)", data_dir, xls_glob, len(xls_files))
    if not xls_files:
        log.error("No matching XLS files found. Aborting.")
        sys.exit(1)

    log.info("\nFound %d metadata file(s):", len(xls_files))
    xls_stats: list = []
    frames:    list = []

    for xls_path in xls_files:
        df = load_metadata(xls_path)
        xls_stats.append({"file": xls_path.name, "rows_after_load": len(df)})
        if not df.empty:
            frames.append(df)

    if not frames:
        log.error("All XLS files were empty or failed to load. Aborting.")
        sys.exit(1)

    total_before = sum(len(f) for f in frames)

    # Step 2: Concatenate + cross-file deduplication -------------------------
    # "first" = alphabetically first filename (sorted() above), then original
    # row order within that file. This contract is explicit and deterministic.
    combined   = pd.concat(frames, ignore_index=True)
    cross_dups = combined.duplicated(subset="Isolate_Id", keep=False)
    if cross_dups.any():
        n_dup        = int(cross_dups.sum())
        n_unique_dup = int(combined.loc[cross_dups, "Isolate_Id"].nunique())
        log.warning(
            "\nCross-file duplicates: %d rows share %d Isolate_Id(s)"
            " -- keeping first (alphabetically first XLS, first row in that file)",
            n_dup, n_unique_dup,
        )
        combined = combined.drop_duplicates(subset="Isolate_Id", keep="first")
    else:
        log.info("\nNo cross-file duplicate Isolate_Id(s)")

    log.info(
        "Combined (deduplicated): %d isolates (from %d across all files)",
        len(combined), total_before,
    )
    after_merge = len(combined)

    # Step 3: Discover + parse all FASTA files --------------------------------
    fasta_glob  = "gisaid_epiflu_sequence_*.fasta"
    fasta_files = sorted(data_dir.glob(fasta_glob))
    log.info("\nFASTA glob     : %s/%s  -> %d file(s)", data_dir, fasta_glob, len(fasta_files))
    xls_tags    = {p.stem.replace("gisaid_epiflu_isolates_", "") for p in xls_files}

    fasta_stats: list = []
    all_seqs:    list = []

    for fasta_path in fasta_files:
        tag = fasta_path.stem.replace("gisaid_epiflu_sequence_", "")
        if tag not in xls_tags:
            log.info(
                "  Note: [%s] has no matching XLS; sequences joined by EPI accession",
                fasta_path.name,
            )
        df_seq = parse_fasta_lengths(fasta_path, seqkit)
        fasta_stats.append({
            "file":      fasta_path.name,
            "sequences": len(df_seq),
            "seg_types": int(df_seq["Segment"].nunique()) if not df_seq.empty else 0,
        })
        if not df_seq.empty:
            all_seqs.append(df_seq)

    # Step 4: Join sequence lengths -------------------------------------------
    if all_seqs:
        fasta_df = pd.concat(all_seqs, ignore_index=True).drop_duplicates("EPI_ID")
        log.info("\nTotal unique sequences across all FASTAs: %d", len(fasta_df))
        combined = attach_sequence_lengths(combined, fasta_df)
    else:
        log.warning("No FASTA files found -- all _length columns will be NaN.")
        for seg in CORE_SEGMENTS:
            combined[f"{seg}_length"] = np.nan

    # Step 5: Completeness flags ----------------------------------------------
    combined = add_sequence_completeness(combined)

    # Step 6: Sort for deterministic output -----------------------------------
    combined = combined.sort_values("Isolate_Id", ignore_index=True)

    # Step 7: Summary ---------------------------------------------------------
    africa = combined[combined["Continent"] == "Africa"]
    mdg    = combined[combined["Madagascar"]]
    print_summary("Madagascar", mdg)
    print_summary("All Africa", africa)

    # Step 8: Write output ----------------------------------------------------
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, sep="\t", index=False)
    log.info("\nSaved: %s  (%d rows, %d columns)", out_path, len(combined), len(combined.columns))

    # Step 9: Provenance ------------------------------------------------------
    write_provenance(
        out_path, seqkit,
        xls_stats=xls_stats,
        fasta_stats=fasta_stats,
        stage_counts={
            "total_loaded":      total_before,
            "after_merge_dedup": after_merge,
            "final":             len(combined),
            "madagascar":        len(mdg),
            "africa":            len(africa),
        },
        min_seq_len=MIN_SEQ_LEN,
    )


if __name__ == "__main__":
    main()
