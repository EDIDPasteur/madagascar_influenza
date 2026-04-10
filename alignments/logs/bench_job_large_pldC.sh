#!/usr/bin/env bash
#SBATCH --job-name=bench_large_mafft
#SBATCH --output=/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H1N1.africa_%j.out
#SBATCH --error=/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H1N1.africa_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=seqbio

source /opt/gensoft/adm/Modules/5.6.1/init/bash
module load fasta ruby mafft/7.526

echo "=== BENCHMARK: large tier ==="
echo "Input:    /pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/split/HA_H1N1.africa.fasta"
echo "Seqs:     8621"
echo "Mode:     --auto --nuc --thread 16"
echo "Alloc:    16 CPUs, 32G RAM, 12:00:00 wall"
echo ""

/usr/bin/time -v   mafft --auto --nuc --thread 16 "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/split/HA_H1N1.africa.fasta" > "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/benchmark/HA_H1N1.africa.bench.aln.fasta"   2> "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H1N1.africa_${SLURM_JOB_ID}.timev"

echo ""
echo "=== /usr/bin/time -v output ==="
cat "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H1N1.africa_${SLURM_JOB_ID}.timev"
