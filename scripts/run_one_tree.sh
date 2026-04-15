#!/usr/bin/env bash
# run_one_tree.sh
#
# Runs IQ-TREE (HKY+G, UFBoot2 -B 1000) on a full (no subsampling) FASTA
# via srun. Called by parallel in run_trees.sh.
#
# Two input types are handled automatically based on the filename:
#
#   <SEG>_<SUBTYPE>.africa.full.fasta   (HA / NA per-subtype)
#       → already aligned; runs IQ-TREE directly.
#
#   <SEG>.africa.full.fasta              (internal segments: PB2 PB1 PA NP MP NS)
#       → unaligned merged file; runs MAFFT --auto first, then IQ-TREE.
#
# Usage: bash scripts/run_one_tree.sh <INPUT_FASTA> <TREE_DIR> <LOG_DIR>

set -euo pipefail

INPUT="$1"
TREE_DIR="$2"
LOG_DIR="$3"

NAME=$(basename "${INPUT}" .africa.full.fasta)
PREFIX="${TREE_DIR}/${NAME}"
SUMMARY="${LOG_DIR}/tree_summary.tsv"

# Detect whether MAFFT re-alignment is required.
# Internal merged segments have names like PB2, PB1, PA, NP, MP, NS (no underscore subtype).
NEEDS_REALIGN=false
if [[ "${NAME}" =~ ^(PB2|PB1|PA|NP|MP|NS)$ ]]; then
    NEEDS_REALIGN=true
fi

# Skip if already done
if [[ -f "${PREFIX}.treefile" ]]; then
    echo "[$(date +%H:%M:%S)] SKIP ${NAME} -- tree exists"
    exit 0
fi

N_SEQS=$(grep -c "^>" "${INPUT}")

# Resource tiers (UFBoot2, no subsampling — trees can be large)
#   ≤ 500 seqs              : 4 CPUs,  4G,  2h
#   501–2000 seqs           : 8 CPUs,  16G, 8h
#   2001–5000 seqs          : 16 CPUs, 32G, 24h
#   > 5000 seqs             : 32 CPUs, 64G, 48h
#   internal (realign)      : 32 CPUs, 128G, 72h  (~13k seqs after MAFFT)
if [[ "${NEEDS_REALIGN}" == "true" ]]; then
    CPUS=32; MEM="128G"; TIME="72:00:00"
elif [[ "${N_SEQS}" -le 500 ]]; then
    CPUS=4;  MEM="4G";   TIME="2:00:00"
elif [[ "${N_SEQS}" -le 2000 ]]; then
    CPUS=8;  MEM="16G";  TIME="8:00:00"
elif [[ "${N_SEQS}" -le 5000 ]]; then
    CPUS=16; MEM="32G";  TIME="24:00:00"
else
    CPUS=32; MEM="64G";  TIME="48:00:00"
fi

echo "[$(date +%H:%M:%S)] START ${NAME} (n=${N_SEQS}, cpus=${CPUS}, mem=${MEM}, realign=${NEEDS_REALIGN})"

if [[ "${NEEDS_REALIGN}" == "true" ]]; then
    # Internal merged segment: MAFFT re-alignment then IQ-TREE with UFBoot2
    ALIGNED="${INPUT%.fasta}.aln.fasta"
    srun --partition=seqbio --cpus-per-task="${CPUS}" --mem="${MEM}" --time="${TIME}" \
        bash -c "
            source /opt/gensoft/adm/Modules/5.6.1/init/bash
            module load fasta ruby mafft/7.526 IQ-TREE/3.1.0
            mafft --auto --nuc --thread ${CPUS} '${INPUT}' > '${ALIGNED}'
            iqtree3 -s '${ALIGNED}' -m HKY+G -B 1000 -T ${CPUS} --seed 42 \
                    --prefix '${PREFIX}' --redo --quiet
        "
else
    # HA / NA per-subtype: already aligned, run IQ-TREE with UFBoot2 directly
    srun --partition=seqbio --cpus-per-task="${CPUS}" --mem="${MEM}" --time="${TIME}" \
        bash -c "
            source /opt/gensoft/adm/Modules/5.6.1/init/bash
            module load IQ-TREE/3.1.0
            iqtree3 -s '${INPUT}' -m HKY+G -B 1000 -T ${CPUS} --seed 42 \
                    --prefix '${PREFIX}' --redo --quiet
        "
fi

STATUS="FAIL"
[[ -f "${PREFIX}.treefile" ]] && STATUS="OK"
printf "%s\t%d\t%s\t%s\n" "${NAME}" "${N_SEQS}" "${NEEDS_REALIGN}" "${STATUS}" >> "${SUMMARY}"
echo "[$(date +%H:%M:%S)] DONE  ${NAME} (${STATUS})"
