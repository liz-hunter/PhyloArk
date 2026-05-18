#!/usr/bin/env bash
#SBATCH -J skani
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --time=24:00:00
#SBATCH --mem=150G

# Assumes skani is pre-compiled and available in path

module load rust/1.91/

skani triangle \
    --ql enterobacteriaceae_clean.txt \
    --ri enterobacteriaceae_clean.txt \
    -t 20 \
    -o skani_triangle.tsv

echo "FINISHED"
date