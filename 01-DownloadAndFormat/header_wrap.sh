#!/usr/bin/env bash
#SBATCH -J fixhead
#SBATCH -c 1
#SBATCH --output=%x_%j.out

module load python/3.12.5 

echo "START"
date

python headers2taxid.py $MYPATH/CMSC701/CMSC701_final/refseq_fastas $MYPATH/CMSC701/CMSC701_final/taxonomy/refseq_accession2taxid.map

echo "END"
date
