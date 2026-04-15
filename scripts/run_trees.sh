#!/usr/bin/env bash
# run_trees.sh
#
# Builds ML phylogenies (IQ-TREE, HKY+G) for trees that contain enough
# Madagascar sequences to support monophyletic-clade analysis:
#
#   Filter  : Mdg >= 10 sequences  AND  non-Mdg Africa >= 20 sequences
#
#   HA / NA : one tree per subtype, from per-subtype alignments.
#             Passing subtypes: H3N2, H1N1, B, H9N2  (4 HA + 4 NA = 8 trees)
#
#   PB2 PB1 PA NP MP NS : one tree per segment, merging ALL subtypes.
#             Internal segments do not cluster by HxNy subtype (Poon 2024).
#             Sequences drawn from alignments/split/ and re-aligned by MAFFT.
#             (6 trees)
#
# Subsampling strategy (all trees):
#   - All Madagascar sequences kept in full
#   - Africa (non-Mdg) subsampled proportionally by country to
#     4 × n_Madagascar  (keeps tree Africa-dominant for clade detection)
#   - Random seed 42 for reproducibility
#
# Output:
#   trees/full/<NAME>.africa.full.fasta          — full inputs (no subsampling)
#   trees/<NAME>.treefile                        — ML trees with UFBoot2 support
#   trees/logs/tree_summary.tsv                  — per-tree status
#
# Usage (from project root):
#   bash scripts/run_trees.sh [--dry-run] [--batch-size N]
#
# --batch-size N : submit at most N new jobs per invocation (default: unlimited).
#                  Re-run to submit the next batch; completed trees are skipped.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALN_DIR="${PROJECT_ROOT}/alignments/aligned"
SPLIT_DIR="${PROJECT_ROOT}/alignments/split"
TREE_DIR="${PROJECT_ROOT}/trees"
FULL_DIR="${TREE_DIR}/full"
LOG_DIR="${TREE_DIR}/logs"

DRY_RUN=false
BATCH_SIZE=0   # 0 = unlimited

ARGS=("$@")
for (( i=0; i<${#ARGS[@]}; i++ )); do
    case "${ARGS[$i]}" in
        --dry-run)    DRY_RUN=true ;;
        --batch-size) BATCH_SIZE="${ARGS[$((i+1))]}"; i=$((i+1)) ;;
    esac
done

source /opt/gensoft/adm/Modules/5.6.1/init/bash
module load parallel/20200222

mkdir -p "${FULL_DIR}" "${LOG_DIR}"

# ---------------------------------------------------------------------------
# 1. Prepare input FASTAs (all sequences — no subsampling)
#    Phase 1A: HA and NA — per-subtype, from aligned files
#    Phase 1B: Internal segments — merged across all subtypes, from split files
# ---------------------------------------------------------------------------
echo "[1/2] Preparing input FASTAs (all sequences, no subsampling)..."

export ALN_DIR SPLIT_DIR FULL_DIR

python3 - <<'PYEOF'
import os, re
from pathlib import Path

ALN_DIR       = Path(os.environ["ALN_DIR"])
SPLIT_DIR     = Path(os.environ["SPLIT_DIR"])
FULL_DIR      = Path(os.environ["FULL_DIR"])
MDG_MIN       = 10   # skip tree if fewer Madagascar sequences
AFR_MIN       = 20   # skip tree if fewer Africa sequences
INTERNAL_SEGS = {"PB2", "PB1", "PA", "NP", "MP", "NS"}

def is_madagascar(header):
    if '|' in header:
        fields = header[1:].split('|')
        return len(fields) > 2 and 'Madagascar' in fields[2]
    return 'Madagascar' in header

def read_fasta(path, strip_gaps=False):
    records, header, seq_lines = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    seq = ''.join(seq_lines)
                    if strip_gaps:
                        seq = seq.replace('-', '')
                    records.append((header, seq))
                header, seq_lines = line, []
            elif header is not None:
                seq_lines.append(line)
    if header is not None:
        seq = ''.join(seq_lines)
        if strip_gaps:
            seq = seq.replace('-', '')
        records.append((header, seq))
    return records

def write_fasta(path, records):
    with open(path, 'w') as fh:
        for header, seq in records:
            fh.write(header + '\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')

written = 0
skipped = 0
print(f"  {'Name':<30} {'Source':<12} {'Total':>7}  {'Mdg':>4}  {'Afr':>4}")
print(f"  {'-'*30} {'-'*12} {'-'*7}  {'-'*4}  {'-'*4}")

# Phase 1A: HA and NA — per-subtype from aligned files (all sequences, gaps kept)
for aln_file in sorted(ALN_DIR.glob('*.africa.aln.fasta')):
    name = re.sub(r'\.africa\.aln$', '', aln_file.stem)
    seg  = name.split('_')[0]
    if seg not in ('HA', 'NA'):
        continue
    records   = read_fasta(aln_file, strip_gaps=False)
    n_mdg_tot = sum(1 for h, _ in records if     is_madagascar(h))
    n_afr_tot = sum(1 for h, _ in records if not is_madagascar(h))
    if n_mdg_tot < MDG_MIN or n_afr_tot < AFR_MIN:
        print(f"  {name:<30} {'SKIP':<12} {len(records):>7,}  {n_mdg_tot:>4}  {n_afr_tot:>4}")
        skipped += 1
        continue
    out = FULL_DIR / f"{name}.africa.full.fasta"
    write_fasta(out, records)
    print(f"  {name:<30} {'aligned':<12} {len(records):>7,}  {n_mdg_tot:>4}  {n_afr_tot:>4}")
    written += 1

# Phase 1B: internal segments — merge ALL subtypes from unaligned split files
# (will be re-aligned by MAFFT inside run_one_tree.sh)
for seg in sorted(INTERNAL_SEGS):
    split_files = sorted(SPLIT_DIR.glob(f'{seg}_*.africa.fasta'))
    if not split_files:
        continue
    all_records = []
    for sf in split_files:
        all_records.extend(read_fasta(sf, strip_gaps=False))
    n_mdg_tot = sum(1 for h, _ in all_records if     is_madagascar(h))
    n_afr_tot = sum(1 for h, _ in all_records if not is_madagascar(h))
    if n_mdg_tot < MDG_MIN or n_afr_tot < AFR_MIN:
        print(f"  {seg:<30} {'SKIP':<12} {len(all_records):>7,}  {n_mdg_tot:>4}  {n_afr_tot:>4}")
        skipped += 1
        continue
    out = FULL_DIR / f"{seg}.africa.full.fasta"   # no subtype → signals realign needed
    write_fasta(out, all_records)
    print(f"  {seg:<30} {'split+merge':<12} {len(all_records):>7,}  {n_mdg_tot:>4}  {n_afr_tot:>4}")
    written += 1

print(f"\n  Written: {written}  Skipped (filter): {skipped}")
print(f"  Output: {FULL_DIR}/")
PYEOF

# ---------------------------------------------------------------------------
# 2. Run IQ-TREE via parallel + srun (one srun per tree, -j N concurrent)
# ---------------------------------------------------------------------------
echo ""
echo "[2/2] Running IQ-TREE jobs (UFBoot2 -B 1000) via parallel + srun..."

JOBS="${BATCH_SIZE}"
[[ "${BATCH_SIZE}" -eq 0 ]] && JOBS=16   # unlimited → all at once (16 trees)

# Build list of pending (not yet tree'd) FASTA files
PENDING=()
for INPUT in "${FULL_DIR}"/*.africa.full.fasta; do
    NAME=$(basename "${INPUT}" .africa.full.fasta)
    if [[ -f "${TREE_DIR}/${NAME}.treefile" ]]; then
        continue
    fi
    PENDING+=("${INPUT}")
done
TOTAL="${#PENDING[@]}"

if [[ "${TOTAL}" -eq 0 ]]; then
    echo "All trees already done."
    exit 0
fi

echo "  ${TOTAL} pending, ${JOBS} concurrent srun slots."

if "${DRY_RUN}"; then
    echo "  [dry-run] Would submit ${TOTAL} IQ-TREE jobs."
    exit 0
fi

printf '%s\n' "${PENDING[@]}" | \
    parallel -j "${JOBS}" \
        --joblog "${LOG_DIR}/parallel_trees.log" \
        --resume \
        "bash ${PROJECT_ROOT}/scripts/run_one_tree.sh {} ${TREE_DIR} ${LOG_DIR}"

echo ""
echo "Done. Trees in ${TREE_DIR}/"
echo "Summary: ${LOG_DIR}/tree_summary.tsv"
