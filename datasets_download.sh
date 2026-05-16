#!/usr/bin/env bash
#SBATCH -J datasets
#SBATCH -N 1
#SBATCH -c 1

echo "START"
date

module load ncbi-datasets-cli/16.27

cd $MYPATH/CMSC701/CMSC701_final

sleep $(( (RANDOM % 5 + 1) * 60 ))

echo "BEGIN DOWNLOAD"
datasets download genome accession --dehydrated --inputfile refseq_accessions.txt --filename refseq.zip

echo "BEGIN UNZIP"
unzip refseq.zip -d refseq

echo "REHYDRATE"
datasets rehydrate --directory refseq

echo "CHECK MD5"
cd refseq
md5sum -c md5sum.txt > ../checks_refseq.out

echo "END"
date