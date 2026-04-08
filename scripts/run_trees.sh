#!/usr/bin/env bash
# run_trees.sh
#
# For each africa-scoped alignment, build an ML phylogeny (IQ-TREE, HKY+G)
# on a subsampled dataset:
#   - All Madagascar sequences (kept in full)
#   - Up to AFRICA_TARGET African sequences (proportional by country, seed 42)
#
# Output:
#   trees/subsampled/<SEG>_<SUBTYPE>.africa.subsampled.fasta  — inputs
#   trees/<SEG>_<SUBTYPE>.treefile                            — ML trees
#   trees/logs/tree_summary.tsv                               — per-tree status
#
# Usage (from project root):
#   bash scripts/run_trees.sh [--dry-run] [--batch-size N]
#
# --batch-size N : submit at most N new jobs per invocation (default: unlimited).
#                  Re-run to submit the next batch; completed trees are skipped.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALN_DIR="${PROJECT_ROOT}/alignments/aligned"
TREE_DIR="${PROJECT_ROOT}/trees"
SUB_DIR="${TREE_DIR}/subsampled"
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

mkdir -p "${SUB_DIR}" "${LOG_DIR}"

# ---------------------------------------------------------------------------
# 1. Subsample each africa alignment: all Madagascar + proportional Africa
# ---------------------------------------------------------------------------
echo "[1/2] Subsampling alignments for tree-building..."

export ALN_DIR SUB_DIR

python3 - <<'PYEOF'
import os, re, random
from pathlib import Path
from collections import defaultdict

random.seed(42)

ALN_DIR       = Path(os.environ["ALN_DIR"])
SUB_DIR       = Path(os.environ["SUB_DIR"])
AFRICA_TARGET = 200   # total non-Madagascar Africa sequences to sample per tree

def is_madagascar(header):
    """Return True if the sequence originates from Madagascar."""
    # GISAID pipe-delimited: >EPI...|SEG|A/country/isolate/year|...
    if '|' in header:
        fields = header[1:].split('|')
        return len(fields) > 2 and 'Madagascar' in fields[2]
    # Norosoa: >A/duck/Madagascar/74316-22/2022_H1N2_PB2
    return 'Madagascar' in header

def get_country(header):
    """Extract country string from a GISAID isolate header."""
    if '|' in header:
        fields = header[1:].split('|')
        if len(fields) > 2:
            parts = fields[2].split('/')
            return parts[1] if len(parts) > 1 else 'Unknown'
    parts = header[1:].split('/')
    return parts[1] if len(parts) > 1 else 'Unknown'

def read_fasta(path):
    """Read (possibly aligned) FASTA; return list of (header, seq) tuples."""
    records, header, seq_lines = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_lines)))
                header, seq_lines = line, []
            elif header is not None:
                seq_lines.append(line)
    if header is not None:
        records.append((header, ''.join(seq_lines)))
    return records

written = 0
print(f"  {'Name':<30} {'Total':>6}  {'Mdg':>4}  {'Afr_sampled':>11}  {'Tree_n':>6}")
print(f"  {'-'*30} {'-'*6}  {'-'*4}  {'-'*11}  {'-'*6}")

for aln_file in sorted(ALN_DIR.glob('*.africa.aln.fasta')):
    name = re.sub(r'\.africa\.aln$', '', aln_file.stem)   # e.g. HA_H3N2
    out  = SUB_DIR / f"{name}.africa.subsampled.fasta"

    records  = read_fasta(aln_file)
    mdg_seqs = [(h, s) for h, s in records if     is_madagascar(h)]
    afr_seqs = [(h, s) for h, s in records if not is_madagascar(h)]

    # Proportional sampling of non-Madagascar Africa by country
    if len(afr_seqs) <= AFRICA_TARGET:
        sampled_afr = afr_seqs
    else:
        by_country = defaultdict(list)
        for h, s in afr_seqs:
            by_country[get_country(h)].append((h, s))
        total_afr   = len(afr_seqs)
        sampled_afr = []
        for country, seqs in by_country.items():
            n_take = max(1, round(len(seqs) / total_afr * AFRICA_TARGET))
            sampled_afr.extend(random.sample(seqs, min(n_take, len(seqs))))

    all_seqs = mdg_seqs + sampled_afr

    with open(out, 'w') as fh:
        for header, seq in all_seqs:
            fh.write(header + '\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')

    print(f"  {name:<30} {len(records):>6}  {len(mdg_seqs):>4}  {len(sampled_afr):>11}  {len(all_seqs):>6}")
    written += 1

print(f"\n  Done -- {written} subsampled FASTAs written to {SUB_DIR}/")
PYEOF

# ---------------------------------------------------------------------------
# 2. Run IQ-TREE via parallel + srun (one srun per tree, -j N concurrent)
# ---------------------------------------------------------------------------
echo ""
echo "[2/2] Running IQ-TREE jobs via parallel + srun..."

JOBS="${BATCH_SIZE}"
[[ "${BATCH_SIZE}" -eq 0 ]] && JOBS=265   # unlimited → all at once

# Build list of pending (not yet tree'd) FASTA files
PENDING=()
for INPUT in "${SUB_DIR}"/*.africa.subsampled.fasta; do
    NAME=$(basename "${INPUT}" .africa.subsampled.fasta)
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
