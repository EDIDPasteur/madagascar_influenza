#!/usr/bin/env bash
# parse_benchmark.sh
#
# After all 4 benchmark jobs finish, parse /usr/bin/time -v output and
# sacct stats to print a calibration table and suggest revised tier thresholds.
#
# Usage (from project root):
#   bash scripts/parse_benchmark.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${PROJECT_ROOT}/alignments/logs"
JOB_ID_FILE="${LOG_DIR}/benchmark_jobs.txt"

if [[ ! -f "${JOB_ID_FILE}" ]]; then
    echo "ERROR: ${JOB_ID_FILE} not found. Run benchmark_alignment.sh first."
    exit 1
fi

echo ""
echo "=== Benchmark results ==="
echo ""
printf "%-8s %-30s %6s %5s %8s %10s %10s %6s\n" \
    "TIER" "FILE" "N_SEQS" "CPUs" "MEM_REQ" "PEAK_RAM" "WALL_TIME" "JOB_ID"
printf "%-8s %-30s %6s %5s %8s %10s %10s %6s\n" \
    "--------" "------------------------------" "------" "-----" "--------" "----------" "----------" "------"

while IFS=: read -r JOB_ID TIER NAME N_SEQS CPUS MEM; do
    # Find .timev file for this job
    TIMEV_FILE=$(ls "${LOG_DIR}/bench_${NAME}_${JOB_ID}.timev" 2>/dev/null || echo "")

    PEAK_RAM="N/A"
    WALL_TIME="N/A"

    if [[ -n "${TIMEV_FILE}" && -f "${TIMEV_FILE}" ]]; then
        # Parse /usr/bin/time -v output
        PEAK_KB=$(grep "Maximum resident set size" "${TIMEV_FILE}" | awk '{print $NF}')
        if [[ -n "${PEAK_KB}" ]]; then
            PEAK_MB=$(( PEAK_KB / 1024 ))
            PEAK_GB=$(echo "scale=1; ${PEAK_MB} / 1024" | bc)
            PEAK_RAM="${PEAK_GB}G"
        fi
        WALL_RAW=$(grep "Elapsed (wall clock) time" "${TIMEV_FILE}" | awk '{print $NF}')
        [[ -n "${WALL_RAW}" ]] && WALL_TIME="${WALL_RAW}"
    fi

    # Fallback: sacct
    if [[ "${PEAK_RAM}" == "N/A" ]] && command -v sacct &>/dev/null; then
        SACCT_MEM=$(sacct -j "${JOB_ID}" --format=MaxRSS --noheader 2>/dev/null | \
                    grep -v "^$\|\.bat" | head -1 | tr -d ' ')
        [[ -n "${SACCT_MEM}" ]] && PEAK_RAM="${SACCT_MEM}"
        SACCT_TIME=$(sacct -j "${JOB_ID}" --format=Elapsed --noheader 2>/dev/null | \
                     grep -v "^$\|\.bat" | head -1 | tr -d ' ')
        [[ -n "${SACCT_TIME}" ]] && WALL_TIME="${SACCT_TIME}"
    fi

    printf "%-8s %-30s %6s %5s %8s %10s %10s %6s\n" \
        "${TIER}" "${NAME}" "${N_SEQS}" "${CPUS}" "${MEM}" "${PEAK_RAM}" "${WALL_TIME}" "${JOB_ID}"

done < "${JOB_ID_FILE}"

echo ""
echo "Tip: if PEAK_RAM > MEM_REQ, increase that tier's memory in align_segments.sh"
echo "Tip: if WALL_TIME approaches the limit, increase --time for that tier"
echo "Tip: if PEAK_RAM << MEM_REQ/2, you can reduce memory (better queue priority)"
echo ""
echo "To check job status: squeue -u \$USER"
echo "To check completed:  sacct -u \$USER --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem -S today"
