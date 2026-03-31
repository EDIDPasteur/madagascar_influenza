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

# Print reportseff summary for all job IDs
JOB_IDS_CSV=$(awk -F: '{printf "%s,",$1}' "${JOB_ID_FILE}" | sed 's/,$//')
reportseff ${JOB_IDS_CSV} 2>/dev/null || \
    echo "reportseff not available — falling back to sacct:"

echo ""
echo "=== Per-tier detail ==="
echo ""
printf "%-8s %-30s %6s %5s %8s %10s %10s %6s\n" \
    "TIER" "FILE" "N_SEQS" "CPUs" "MEM_REQ" "PEAK_RAM" "WALL_TIME" "JOB_ID"
printf "%-8s %-30s %6s %5s %8s %10s %10s %6s\n" \
    "--------" "------------------------------" "------" "-----" "--------" "----------" "----------" "------"

while IFS=: read -r JOB_ID TIER NAME N_SEQS CPUS MEM; do
    PEAK_RAM="N/A"
    WALL_TIME="N/A"
    EXIT_STATUS="N/A"

    # /usr/bin/time -v output embedded in .out log
    OUT_FILE=$(ls "${LOG_DIR}/bench_${NAME}_${JOB_ID}.out" 2>/dev/null || echo "")
    if [[ -n "${OUT_FILE}" && -f "${OUT_FILE}" ]]; then
        PEAK_KB=$(grep "Maximum resident set size" "${OUT_FILE}" | awk '{print $NF}')
        if [[ -n "${PEAK_KB}" ]]; then
            PEAK_GB=$(echo "scale=2; ${PEAK_KB} / 1024 / 1024" | bc)
            PEAK_RAM="${PEAK_GB}G"
        fi
        WALL_RAW=$(grep "Elapsed (wall clock) time" "${OUT_FILE}" | awk '{print $NF}')
        [[ -n "${WALL_RAW}" ]] && WALL_TIME="${WALL_RAW}"
        EXIT_RAW=$(grep "Exit status" "${OUT_FILE}" | awk '{print $NF}')
        [[ -n "${EXIT_RAW}" ]] && EXIT_STATUS="${EXIT_RAW}"
    fi

    STATUS_FLAG=""
    [[ "${EXIT_STATUS}" != "0" && "${EXIT_STATUS}" != "N/A" ]] && STATUS_FLAG=" !! EXIT=${EXIT_STATUS}"

    printf "%-8s %-30s %6s %5s %8s %10s %10s %6s%s\n" \
        "${TIER}" "${NAME}" "${N_SEQS}" "${CPUS}" "${MEM}" \
        "${PEAK_RAM}" "${WALL_TIME}" "${JOB_ID}" "${STATUS_FLAG}"

done < "${JOB_ID_FILE}"

echo ""
echo "Tip: if PEAK_RAM > MEM_REQ, increase that tier's memory in align_segments.sh"
echo "Tip: if WALL_TIME approaches the limit, increase --time for that tier"
echo "Tip: if PEAK_RAM << MEM_REQ/2, you can safely reduce memory (better queue priority)"
echo "Tip: jobs marked !! EXIT≠0 failed — check log: alignments/logs/bench_<NAME>_<JOBID>.out"
