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
#   trees/subsampled/<NAME>.africa.subsampled.fasta  — subsampled inputs
#   trees/<NAME>.treefile                            — ML trees
#   trees/logs/tree_summary.tsv                      — per-tree status
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
# 1. Subsample alignments for tree-building
#    Phase 1A: HA and NA — per-subtype, from aligned files
#    Phase 1B: Internal segments — merged across all subtypes, from split files
# ---------------------------------------------------------------------------
echo "[1/2] Subsampling alignments for tree-building..."

export ALN_DIR SPLIT_DIR SUB_DIR

python3 - <<'PYEOF'
import os, re, random
from pathlib import Path
from collections import defaultdict

random.seed(42)

ALN_DIR           = Path(os.environ["ALN_DIR"])
SPLIT_DIR         = Path(os.environ["SPLIT_DIR"])
SUB_DIR           = Path(os.environ["SUB_DIR"])
AFRICA_MULTIPLIER = 4    # Africa target = 4 × n_Madagascar per tree
MDG_MIN           = 10   # skip tree if fewer Madagascar sequences
AFR_MIN           = 20   # skip tree if fewer Africa sequences
INTERNAL_SEGS     = {"PB2", "PB1", "PA", "NP", "MP", "NS"}

def is_madagascar(header):
    if '|' in header:
        fields = header[1:].split('|')
        return len(fields) > 2 and 'Madagascar' in fields[2]
    return 'Madagascar' in header

def get_country(header):
    if '|' in header:
        fields = header[1:].split('|')
        if len(fields) > 2:
            parts = fields[2].split('/')
            return parts[1] if len(parts) > 1 else 'Unknown'
    return 'Unknown'

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

def subsample_africa(records):
    """Keep all Madagascar; proportionally subsample Africa to 4 × n_Madagascar."""
    mdg_seqs = [(h, s) for h, s in records if     is_madagascar(h)]
    afr_seqs = [(h, s) for h, s in records if not is_madagascar(h)]
    africa_target = len(mdg_seqs) * AFRICA_MULTIPLIER
    if len(afr_seqs) <= africa_target:
        return mdg_seqs + afr_seqs, len(mdg_seqs), len(afr_seqs)
    by_country = defaultdict(list)
    for h, s in afr_seqs:
        by_country[get_country(h)].append((h, s))
    total_afr = len(afr_seqs)
    sampled_afr = []
    for country, seqs in by_country.items():
        n_take = max(1, round(len(seqs) / total_afr * africa_target))
        sampled_afr.extend(random.sample(seqs, min(n_take, len(seqs))))
    return mdg_seqs + sampled_afr, len(mdg_seqs), len(sampled_afr)

def write_fasta(path, records):
    with open(path, 'w') as fh:
        for header, seq in records:
            fh.write(header + '\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')

written = 0
print(f"  {'Name':<30} {'Source':<10} {'Total':>7}  {'Mdg':>4}  {'Afr':>4}  {'Tree_n':>6}")
print(f"  {'-'*30} {'-'*10} {'-'*7}  {'-'*4}  {'-'*4}  {'-'*6}")

skipped = 0

# Phase 1A: HA and NA — per-subtype from aligned files (gaps kept, already aligned)
for aln_file in sorted(ALN_DIR.glob('*.africa.aln.fasta')):
    name = re.sub(r'\.africa\.aln$', '', aln_file.stem)
    seg  = name.split('_')[0]
    if seg not in ('HA', 'NA'):
        continue
    records    = read_fasta(aln_file, strip_gaps=False)
    n_mdg_tot  = sum(1 for h, _ in records if     is_madagascar(h))
    n_afr_tot  = sum(1 for h, _ in records if not is_madagascar(h))
    if n_mdg_tot < MDG_MIN or n_afr_tot < AFR_MIN:
        print(f"  {name:<30} {'SKIP':<11} {len(records):>7,}  {n_mdg_tot:>4}  {n_afr_tot:>4}  {'---':>6}")
        skipped += 1
        continue
    out = SUB_DIR / f"{name}.africa.subsampled.fasta"
    all_seqs, n_mdg, n_afr = subsample_africa(records)
    write_fasta(out, all_seqs)
    print(f"  {name:<30} {'aligned':<11} {len(records):>7,}  {n_mdg:>4}  {n_afr:>4}  {len(all_seqs):>6}")
    written += 1

# Phase 1B: internal segments — merge ALL subtypes from unaligned split files
# These will be re-aligned by MAFFT inside run_one_tree.sh
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
        print(f"  {seg:<30} {'SKIP':<11} {len(all_records):>7,}  {n_mdg_tot:>4}  {n_afr_tot:>4}  {'---':>6}")
        skipped += 1
        continue
    out = SUB_DIR / f"{seg}.africa.subsampled.fasta"  # no subtype → signals realign needed
    all_seqs, n_mdg, n_afr = subsample_africa(all_records)
    write_fasta(out, all_seqs)
    print(f"  {seg:<30} {'split+merge':<11} {len(all_records):>7,}  {n_mdg:>4}  {n_afr:>4}  {len(all_seqs):>6}")
    written += 1

print(f"\n  Written: {written}  Skipped (filter): {skipped}")
print(f"  Output: {SUB_DIR}/")
PYEOF

# ---------------------------------------------------------------------------
# 2. Run IQ-TREE via parallel + srun (one srun per tree, -j N concurrent)
# ---------------------------------------------------------------------------
echo ""
echo "[2/2] Running IQ-TREE jobs via parallel + srun..."

JOBS="${BATCH_SIZE}"
[[ "${BATCH_SIZE}" -eq 0 ]] && JOBS=16   # unlimited → all at once (16 trees)

# Build list of pending (not yet tree'd) FASTA files
# Both *.africa.subsampled.fasta (HA/NA) and *.africa.subsampled.fasta (internal)
PENDING=()
for INPUT in "${SUB_DIR}"/*.africa.subsampled.fasta; do
    # Derive tree name: strip .africa.subsampled.fasta suffix
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
