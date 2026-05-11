#!/usr/bin/env bash
#SBATCH -J datasets
#SBATCH -c 1
#SBATCH --output=%x_%A.out

echo "START"
date

module load ncbi-datasets-cli/16.27

cd $MYPATH/CMSC701/CMSC701_final

sleep $(( (RANDOM % 5 + 1) * 60 ))

echo "BEGIN DOWNLOAD"
datasets download genome accession --dehydrated --inputfile enterobacter_genbank.tsv --filename entero_accessions.zip

echo "BEGIN UNZIP"
unzip entero_accessions.zip -d entero

echo "REHYDRATE"
datasets rehydrate --directory entero

echo "CHECK MD5"
cd entero
md5sum -c md5sum.txt > check_entero.out

echo "END"
date
