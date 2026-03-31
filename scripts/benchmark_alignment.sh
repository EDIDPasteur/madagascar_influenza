#!/usr/bin/env bash
# benchmark_alignment.sh
#
# Submit one MAFFT alignment job per resource tier to measure real wall-time
# and peak memory BEFORE launching the full 425-job batch.
#
# Tier targets (all HA, worst-case for memory):
#   tiny   n ~  10 : HA_H11N2.madagascar  (L-INS-i, 4 CPUs, 8G, 4h)
#   small  n ~ 651 : HA_H5N8.africa       (--auto,   4 CPUs, 8G, 4h)
#   medium n ~2303 : HA_H5N1.africa       (--auto,   8 CPUs, 16G, 8h)
#   large  n ~8621 : HA_H1N1.africa       (--auto,  16 CPUs, 32G, 12h)
#
# After all 4 jobs finish, run:
#   bash scripts/parse_benchmark.sh
#
# Usage (from project root):
#   bash scripts/benchmark_alignment.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SPLIT_DIR="${PROJECT_ROOT}/alignments/split"
BENCH_DIR="${PROJECT_ROOT}/alignments/benchmark"
LOG_DIR="${PROJECT_ROOT}/alignments/logs"

mkdir -p "${BENCH_DIR}" "${LOG_DIR}"

# Verify split files exist
N_SPLIT=$(find "${SPLIT_DIR}" -name "*.fasta" 2>/dev/null | wc -l)
if [[ "${N_SPLIT}" -lt 10 ]]; then
    echo "ERROR: split files not found in ${SPLIT_DIR}/"
    echo "       Run: bash scripts/align_segments.sh --dry-run   (this populates alignments/split/)"
    exit 1
fi
echo "Found ${N_SPLIT} split FASTA files in alignments/split/"

# ---------------------------------------------------------------------------
# Tier definitions: name, MAFFT_ARGS, CPUS, MEM, TIME
# ---------------------------------------------------------------------------
declare -A TIER_MAFFT=( [tiny]="--localpair --maxiterate 1000 --nuc --thread 4"
                        [small]="--auto --nuc --thread 4"
                        [medium]="--auto --nuc --thread 8"
                        [large]="--auto --nuc --thread 16" )
declare -A TIER_CPUS=(  [tiny]=4 [small]=4 [medium]=8  [large]=16 )
declare -A TIER_MEM=(   [tiny]="8G" [small]="8G" [medium]="16G" [large]="32G" )
declare -A TIER_TIME=(  [tiny]="04:00:00" [small]="04:00:00" [medium]="08:00:00" [large]="12:00:00" )
declare -A TIER_FILE=(  [tiny]="HA_H11N2.madagascar"
                        [small]="HA_H5N8.africa"
                        [medium]="HA_H5N1.africa"
                        [large]="HA_H1N1.africa" )

# ---------------------------------------------------------------------------
# Submit one job per tier
# ---------------------------------------------------------------------------
JOB_IDS=()
echo ""
echo "Submitting 4 benchmark jobs (partition: seqbio)..."
echo ""

for TIER in tiny small medium large; do
    NAME="${TIER_FILE[$TIER]}"
    INPUT="${SPLIT_DIR}/${NAME}.fasta"

    if [[ ! -f "${INPUT}" ]]; then
        echo "  [WARN] ${TIER}: ${INPUT} not found -- skipping"
        continue
    fi

    N_SEQS=$(grep -c "^>" "${INPUT}")
    OUTPUT="${BENCH_DIR}/${NAME}.bench.aln.fasta"
    LOGFILE="${LOG_DIR}/bench_${NAME}"
    MAFFT_ARGS="${TIER_MAFFT[$TIER]}"
    CPUS="${TIER_CPUS[$TIER]}"
    MEM="${TIER_MEM[$TIER]}"
    TIME="${TIER_TIME[$TIER]}"

    JOB_SCRIPT=$(mktemp "${LOG_DIR}/bench_job_${TIER}_XXXX.sh")
    cat > "${JOB_SCRIPT}" <<SLURM
#!/usr/bin/env bash
#SBATCH --job-name=bench_${TIER}_mafft
#SBATCH --output=${LOGFILE}_%j.out
#SBATCH --error=${LOGFILE}_%j.err
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --partition=seqbio

source /opt/gensoft/adm/Modules/5.6.1/init/bash
module load fasta ruby mafft/7.526

echo "=== BENCHMARK: ${TIER} tier ==="
echo "Input:    ${INPUT}"
echo "Seqs:     ${N_SEQS}"
echo "Mode:     ${MAFFT_ARGS}"
echo "Alloc:    ${CPUS} CPUs, ${MEM} RAM, ${TIME} wall"
echo ""

/usr/bin/time -v \
  mafft ${MAFFT_ARGS} "${INPUT}" > "${OUTPUT}" \
  2> "${LOGFILE}_\${SLURM_JOB_ID}.timev"

echo ""
echo "=== /usr/bin/time -v output ==="
cat "${LOGFILE}_\${SLURM_JOB_ID}.timev"
SLURM

    JOB_ID=$(sbatch "${JOB_SCRIPT}" | awk '{print $NF}')
    JOB_IDS+=("${JOB_ID}:${TIER}:${NAME}:${N_SEQS}:${CPUS}:${MEM}")
    echo "  [${TIER}] n=${N_SEQS}, ${CPUS} CPUs ${MEM} ${TIME} -> job ${JOB_ID}  (${NAME})"
done

# ---------------------------------------------------------------------------
# Write job ID file for parse_benchmark.sh
# ---------------------------------------------------------------------------
JOB_ID_FILE="${LOG_DIR}/benchmark_jobs.txt"
printf "%s\n" "${JOB_IDS[@]}" > "${JOB_ID_FILE}"
echo ""
echo "Job IDs written to: ${JOB_ID_FILE}"

echo ""
echo "Monitor with:  squeue -u \$USER"
echo "When done, run: bash scripts/parse_benchmark.sh"
