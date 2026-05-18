#!/usr/bin/env bash
#SBATCH -J kraken2
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=500G
#SBATCH --output=%x_%A.out

echo "START"
date
module load kraken2/2.1.5

cd $MYPATH/CMSC701/CMSC701_final/

kraken2 --paired --threads 20 --use-names --report-minimizer-data \
--gzip-compressed \
--output kraken2_reports/ani.out \
--report kraken2_reports/ani_report.txt \
--db ani_db \
test_data/nextseq_1M_entero_R1.fastq.gz \
test_data/nextseq_1M_entero_R2.fastq.gz 

echo "END"
date
