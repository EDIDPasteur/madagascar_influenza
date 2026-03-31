#!/usr/bin/env bash
# run_one_alignment.sh
#
# Aligns a single split FASTA file via srun. Called by parallel in align_segments.sh.
#
# Usage: bash scripts/run_one_alignment.sh <INPUT_FASTA> <ALN_DIR> <LOG_DIR>

set -euo pipefail

INPUT="$1"
ALN_DIR="$2"
LOG_DIR="$3"

BASENAME=$(basename "${INPUT}" .fasta)
OUTPUT="${ALN_DIR}/${BASENAME}.aln.fasta"

# Skip if already done (double-check besides parallel --resume)
if [[ -f "${OUTPUT}" ]]; then
    echo "[$(date +%H:%M:%S)] SKIP ${BASENAME} -- already aligned"
    exit 0
fi

N_SEQS=$(grep -c "^>" "${INPUT}")

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

MODE="${MAFFT_ARGS%% *}"
echo "[$(date +%H:%M:%S)] START ${BASENAME} (n=${N_SEQS}, mode=${MODE}, cpus=${CPUS}, mem=${MEM})"

srun --partition=seqbio --cpus-per-task="${CPUS}" --mem="${MEM}" \
    bash -c "
        source /opt/gensoft/adm/Modules/5.6.1/init/bash
        module load fasta ruby mafft/7.526
        mafft ${MAFFT_ARGS} '${INPUT}' > '${OUTPUT}'
    "

# Verify and log
N_OUT=$(grep -c "^>" "${OUTPUT}" 2>/dev/null || echo 0)
SUMMARY="${LOG_DIR}/alignment_summary.tsv"
if [[ "${N_OUT}" -eq "${N_SEQS}" ]]; then
    echo -e "${BASENAME}\t${N_SEQS}\t${N_OUT}\tOK" >> "${SUMMARY}"
    echo "[$(date +%H:%M:%S)] DONE  ${BASENAME} (${N_OUT} seqs)"
else
    echo -e "${BASENAME}\t${N_SEQS}\t${N_OUT}\tWARN_COUNT_MISMATCH" >> "${SUMMARY}"
    echo "[$(date +%H:%M:%S)] WARN  ${BASENAME} input=${N_SEQS} output=${N_OUT}"
fi
