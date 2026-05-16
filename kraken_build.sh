#!/usr/bin/env bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -p tera
#SBATCH --mem=999G

module load kraken2/2.1.5

echo "START"
date

# 1Tb memory available with 50G for overhead
kraken2-build --threads 20 --build --max-db-size 1020054732800 --db $MYPATH/CMSC701/CMSC701_final/taxonomy_db

echo "END"
date
