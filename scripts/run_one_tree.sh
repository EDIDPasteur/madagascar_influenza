#!/usr/bin/env bash
# run_one_tree.sh
#
# Runs IQ-TREE (HKY+G, no branch support) on a single subsampled FASTA
# via srun. Called by parallel in run_trees.sh.
#
# Usage: bash scripts/run_one_tree.sh <INPUT_FASTA> <TREE_DIR> <LOG_DIR>

set -euo pipefail

INPUT="$1"
TREE_DIR="$2"
LOG_DIR="$3"

NAME=$(basename "${INPUT}" .africa.subsampled.fasta)
PREFIX="${TREE_DIR}/${NAME}"
SUMMARY="${LOG_DIR}/tree_summary.tsv"

# Skip if already done
if [[ -f "${PREFIX}.treefile" ]]; then
    echo "[$(date +%H:%M:%S)] SKIP ${NAME} -- tree exists"
    exit 0
fi

N_SEQS=$(grep -c "^>" "${INPUT}")

# Resource tiers for IQ-TREE (subsampled trees are small: max ~575 seqs)
#   ≤ 500 seqs : 4 CPUs, 2G
#   > 500 seqs : 8 CPUs, 4G
if [[ "${N_SEQS}" -le 500 ]]; then
    CPUS=4; MEM="2G"
else
    CPUS=8; MEM="4G"
fi

echo "[$(date +%H:%M:%S)] START ${NAME} (n=${N_SEQS}, cpus=${CPUS}, mem=${MEM})"

srun --partition=seqbio --cpus-per-task="${CPUS}" --mem="${MEM}" \
    bash -c "
        source /opt/gensoft/adm/Modules/5.6.1/init/bash
        module load IQ-TREE/3.1.0
        iqtree3 -s '${INPUT}' -m HKY+G -T ${CPUS} --seed 42 \
                --prefix '${PREFIX}' --redo --quiet
    "

STATUS="FAIL"
[[ -f "${PREFIX}.treefile" ]] && STATUS="OK"
printf "%s\t%d\t%s\n" "${NAME}" "${N_SEQS}" "${STATUS}" >> "${SUMMARY}"
echo "[$(date +%H:%M:%S)] DONE  ${NAME} (${STATUS})"
