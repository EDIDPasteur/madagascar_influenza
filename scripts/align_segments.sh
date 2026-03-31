#!/usr/bin/env bash
# align_segments.sh
#
# For each segment × subtype combination, produce two aligned FASTAs:
#   alignments/aligned/<SEG>_<SUBTYPE>.africa.aln.fasta     (all Africa)
#   alignments/aligned/<SEG>_<SUBTYPE>.madagascar.aln.fasta  (Madagascar only)
#
# MAFFT strategy (applied per file):
#   n ≤ 500  → --localpair --maxiterate 1000  (L-INS-i, most accurate)
#   n > 500  → --auto                         (FFT-NS-2/FFT-NS-1, fast)
#
# Usage (from project root):
#   bash scripts/align_segments.sh [--dry-run] [--batch-size N]
#
# --batch-size N : submit at most N new jobs per invocation (default: unlimited).
#                  Re-run to submit the next batch; completed alignments are skipped.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
SPLIT_DIR="${PROJECT_ROOT}/alignments/split"
ALN_DIR="${PROJECT_ROOT}/alignments/aligned"
LOG_DIR="${PROJECT_ROOT}/alignments/logs"

DRY_RUN=false
BATCH_SIZE=0   # 0 = unlimited

ARGS=("$@")
for (( i=0; i<${#ARGS[@]}; i++ )); do
    case "${ARGS[$i]}" in
        --dry-run)    DRY_RUN=true ;;
        --batch-size) BATCH_SIZE="${ARGS[$((i+1))]}"; i=$((i+1)) ;;
    esac
done

mkdir -p "${SPLIT_DIR}" "${ALN_DIR}" "${LOG_DIR}"

# ---------------------------------------------------------------------------
# 1. Split all FASTA files by segment x subtype x scope (Africa / Madagascar)
# ---------------------------------------------------------------------------
echo "[1/2] Splitting sequences by segment x subtype..."

export DATA_DIR SPLIT_DIR

python3 - <<'PYEOF'
import os, re
from pathlib import Path
from collections import defaultdict

data_dir  = Path(os.environ["DATA_DIR"])
split_dir = Path(os.environ["SPLIT_DIR"])

SEGMENTS = {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"}
MIN_SEQS  = 3   # skip combos with fewer sequences

def norm_subtype(raw):
    """'A_/_H3N2' -> 'H3N2',  'B' -> 'B',  'A_/_H5N1' -> 'H5N1'"""
    s = raw.strip()
    s = re.sub(r'^[AB]_/_', '', s)
    s = re.sub(r'[^A-Za-z0-9]', '_', s)
    return s if s else "Unknown"

def parse_gisaid_header(line):
    """Pipe-delimited GISAID EpiFlu header.
    e.g. >A/duck/Madagascar/78609/2023|HA|Madagascar|2023-07-01|H9N2|...
    Returns (seg, subtype, is_madagascar).
    """
    parts = line[1:].split("|")
    seg     = parts[1].strip() if len(parts) > 1 else None
    subtype = norm_subtype(parts[4]) if len(parts) > 4 else "Unknown"
    is_mdg  = len(parts) > 2 and "Madagascar" in parts[2]
    return seg, subtype, is_mdg

def parse_norosoa_header(line):
    """Norosoa unpublished header format.
    e.g. >A/duck/Madagascar/74316-22/2022_H1N2_PB2
    Last slash-field is YEAR_SUBTYPE_SEGMENT.
    All sequences are from Madagascar -> is_madagascar always True.
    """
    last_field = line[1:].split("/")[-1]   # e.g. "2022_H1N2_PB2"
    parts = last_field.split("_")          # ['2022', 'H1N2', 'PB2']
    seg     = parts[2].strip() if len(parts) >= 3 and parts[2].strip() in SEGMENTS else None
    subtype = norm_subtype(parts[1]) if len(parts) >= 2 else "Unknown"
    return seg, subtype, True

records = defaultdict(lambda: {"africa": [], "madagascar": []})

def process_fasta(fpath, header_parser):
    header = seq_lines = seg = subtype = is_mdg = None
    seq_lines = []

    def flush():
        if header and seg and subtype and seg in SEGMENTS and subtype != "Unknown":
            entry = (header, "".join(seq_lines))
            records[(seg, subtype)]["africa"].append(entry)
            if is_mdg:
                records[(seg, subtype)]["madagascar"].append(entry)

    with open(fpath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                seq_lines = []
                seg, subtype, is_mdg = header_parser(line)
                header = line
            else:
                seq_lines.append(line)
    flush()

# GISAID files
gisaid_files = sorted(data_dir.glob("gisaid_epiflu_sequence_*.fasta"))
if not gisaid_files:
    print("WARNING: no gisaid_epiflu_sequence_*.fasta found in data/")
for fpath in gisaid_files:
    print(f"  Reading {fpath.name} ...", flush=True)
    process_fasta(fpath, parse_gisaid_header)

# Norosoa unpublished avian sequences (all Madagascar)
norosoa_files = sorted(data_dir.glob("norosoa_avian_*.fasta"))
for fpath in norosoa_files:
    print(f"  Reading {fpath.name} (Norosoa, unpublished) ...", flush=True)
    process_fasta(fpath, parse_norosoa_header)

IUPAC = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv-")

def sanitize(seq):
    """Replace any non-IUPAC nucleotide character with N."""
    return "".join(c if c in IUPAC else "N" for c in seq)

written = 0
skipped = 0
sanitized_total = 0
for (seg, subtype), scopes in sorted(records.items()):
    for scope, entries in scopes.items():
        if len(entries) < MIN_SEQS:
            skipped += 1
            continue
        fname = split_dir / f"{seg}_{subtype}.{scope}.fasta"
        with open(fname, "w") as fh:
            for hdr, seq in entries:
                clean = sanitize(seq)
                n_bad = sum(1 for a, b in zip(seq, clean) if a != b)
                if n_bad:
                    sanitized_total += n_bad
                fh.write(hdr + "\n")
                for i in range(0, len(clean), 80):
                    fh.write(clean[i:i+80] + "\n")
        written += 1

if sanitized_total:
    print(f"\n  WARNING: {sanitized_total} non-IUPAC characters replaced with N across all sequences")

print(f"\n  Done -- {written} files written, {skipped} combos skipped (< {MIN_SEQS} seqs)")
print(f"\n  {'Segment':<6} {'Subtype':<12} {'Africa':>8} {'Madagascar':>12}")
print(f"  {'-'*6} {'-'*12} {'-'*8} {'-'*12}")
for (seg, subtype), scopes in sorted(records.items()):
    n_af  = len(scopes["africa"])
    n_mdg = len(scopes["madagascar"])
    if n_af >= MIN_SEQS:
        mdg_str = str(n_mdg) if n_mdg >= MIN_SEQS else f"{n_mdg} (skip)"
        print(f"  {seg:<6} {subtype:<12} {n_af:>8,} {mdg_str:>12}")
PYEOF

# ---------------------------------------------------------------------------
# 2. Submit one SLURM job per FASTA file to align
# ---------------------------------------------------------------------------
echo ""
echo "[2/2] Submitting MAFFT alignment jobs..."
[[ "${BATCH_SIZE}" -gt 0 ]] && echo "      Batch mode: submitting at most ${BATCH_SIZE} jobs."

submitted=0
skipped=0

for INPUT in "${SPLIT_DIR}"/*.fasta; do
    # Stop if batch limit reached
    if [[ "${BATCH_SIZE}" -gt 0 && "${submitted}" -ge "${BATCH_SIZE}" ]]; then
        break
    fi

    BASENAME=$(basename "${INPUT}" .fasta)
    OUTPUT="${ALN_DIR}/${BASENAME}.aln.fasta"

    if [[ -f "${OUTPUT}" ]]; then
        echo "  [SKIP] ${BASENAME} -- already aligned"
        ((skipped++)) || true
        continue
    fi

    N_SEQS=$(grep -c "^>" "${INPUT}" 2>/dev/null || echo 0)

    # Resource tiers — calibrated from benchmarks (HA segment, worst-case memory):
    #   tiny   n <=   500 : L-INS-i, peak 168 MB → 2G (12x)
    #   small  n <=  2000 : --auto,  peak 186 MB → 2G (11x)
    #   medium n <=  8000 : --auto,  peak 468 MB → 4G ( 9x)
    #   large  n >   8000 : --auto,  peak 5.1 GB → 16G (3x)
    if [[ "${N_SEQS}" -le 500 ]]; then
        MAFFT_ARGS="--localpair --maxiterate 1000 --nuc --thread 4"
        CPUS=4; MEM="2G"
    elif [[ "${N_SEQS}" -le 2000 ]]; then
        MAFFT_ARGS="--auto --nuc --thread 4"
        CPUS=4; MEM="2G"
    elif [[ "${N_SEQS}" -le 8000 ]]; then
        MAFFT_ARGS="--auto --nuc --thread 8"
        CPUS=8; MEM="4G"
    else
        MAFFT_ARGS="--auto --nuc --thread 16"
        CPUS=16; MEM="16G"
    fi

    if $DRY_RUN; then
        MODE=$(echo "${MAFFT_ARGS}" | awk '{print $1}')
        echo "  [DRY-RUN] ${BASENAME} (n=${N_SEQS}, mode=${MODE}, cpus=${CPUS}, mem=${MEM})"
        ((submitted++)) || true
        continue
    fi

    JOB_SCRIPT=$(mktemp "${LOG_DIR}/job_${BASENAME}_XXXX.sh")
    cat > "${JOB_SCRIPT}" <<SLURM
#!/usr/bin/env bash
#SBATCH --job-name=mafft_${BASENAME}
#SBATCH --output=${LOG_DIR}/mafft_${BASENAME}_%j.out
#SBATCH --error=${LOG_DIR}/mafft_${BASENAME}_%j.err
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --partition=seqbio

source /opt/gensoft/adm/Modules/5.6.1/init/bash
module load fasta ruby mafft/7.526

echo "Aligning ${BASENAME} (${N_SEQS} sequences, mode=${MAFFT_ARGS%% *})..."
mafft ${MAFFT_ARGS} "${INPUT}" > "${OUTPUT}"

# Verify output and write to summary log
N_OUT=\$(grep -c "^>" "${OUTPUT}" 2>/dev/null || echo 0)
SUMMARY="${LOG_DIR}/alignment_summary.tsv"
if [[ "\${N_OUT}" -eq "${N_SEQS}" ]]; then
    echo -e "${BASENAME}\t${N_SEQS}\t\${N_OUT}\tOK" >> "\${SUMMARY}"
    echo "Done: ${OUTPUT} (\${N_OUT} sequences aligned)"
else
    echo -e "${BASENAME}\t${N_SEQS}\t\${N_OUT}\tWARN_COUNT_MISMATCH" >> "\${SUMMARY}"
    echo "WARNING: input had ${N_SEQS} seqs but output has \${N_OUT}"
fi
SLURM

    JOB_ID=$(sbatch "${JOB_SCRIPT}" | awk '{print $NF}')
    MODE=$(echo "${MAFFT_ARGS}" | awk '{print $1}')
    echo "  Submitted ${BASENAME} (n=${N_SEQS}, mode=${MODE}) -> job ${JOB_ID}"
    ((submitted++)) || true
done

echo ""
if $DRY_RUN; then
    echo "Dry run complete -- ${submitted} jobs would be submitted."
else
    REMAINING=$(( $(ls "${SPLIT_DIR}"/*.fasta | wc -l) - $(ls "${ALN_DIR}"/*.aln.fasta 2>/dev/null | wc -l) - submitted ))
    echo "${submitted} jobs submitted, ${skipped} already done."
    [[ "${BATCH_SIZE}" -gt 0 && "${REMAINING}" -gt 0 ]] && \
        echo "Re-run with --batch-size ${BATCH_SIZE} to submit the next batch (${REMAINING} remaining)."
    echo "Monitor with: squeue -u \$USER"
fi
echo "Alignments will be written to: ${ALN_DIR}/"
