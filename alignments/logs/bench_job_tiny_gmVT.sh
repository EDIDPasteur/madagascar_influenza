#!/usr/bin/env bash
#SBATCH --job-name=bench_tiny_mafft
#SBATCH --output=/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H11N2.madagascar_%j.out
#SBATCH --error=/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H11N2.madagascar_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=seqbio

source /opt/gensoft/adm/Modules/5.6.1/init/bash
module load fasta ruby mafft/7.526

echo "=== BENCHMARK: tiny tier ==="
echo "Input:    /pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/split/HA_H11N2.madagascar.fasta"
echo "Seqs:     10"
echo "Mode:     --localpair --maxiterate 1000 --nuc --thread 4"
echo "Alloc:    4 CPUs, 8G RAM, 04:00:00 wall"
echo ""

/usr/bin/time -v   mafft --localpair --maxiterate 1000 --nuc --thread 4 "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/split/HA_H11N2.madagascar.fasta" > "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/benchmark/HA_H11N2.madagascar.bench.aln.fasta"   2> "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H11N2.madagascar_${SLURM_JOB_ID}.timev"

echo ""
echo "=== /usr/bin/time -v output ==="
cat "/pasteur/appa/scratch/cduitama/EDID/Madagascar_influenza/alignments/logs/bench_HA_H11N2.madagascar_${SLURM_JOB_ID}.timev"
