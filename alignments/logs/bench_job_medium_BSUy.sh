#!/usr/bin/env bash
#SBATCH --job-name=bench_medium_mafft
#SBATCH --output=/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H5N1.africa_%j.out
#SBATCH --error=/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H5N1.africa_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=seqbio

source /opt/gensoft/adm/Modules/5.6.1/init/bash
module load fasta ruby mafft/7.526

echo "=== BENCHMARK: medium tier ==="
echo "Input:    /pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/split/HA_H5N1.africa.fasta"
echo "Seqs:     2303"
echo "Mode:     --auto --nuc --thread 8"
echo "Alloc:    8 CPUs, 16G RAM, 08:00:00 wall"
echo ""

/usr/bin/time -v   mafft --auto --nuc --thread 8 "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/split/HA_H5N1.africa.fasta" > "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/benchmark/HA_H5N1.africa.bench.aln.fasta"   2> "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H5N1.africa_${SLURM_JOB_ID}.timev"

echo ""
echo "=== /usr/bin/time -v output ==="
cat "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H5N1.africa_${SLURM_JOB_ID}.timev"
