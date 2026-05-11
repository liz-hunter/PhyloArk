#!/usr/bin/env bash
#SBATCH -J datasets
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=1-8%2
#SBATCH --output=%x_%a_%A.out

echo "START"
date

module load ncbi-datasets-cli/16.27

cd $MYPATH/CMSC701/CMSC701_final

sleep $(( (RANDOM % 5 + 1) * 60 ))

FILE=$SLURM_ARRAY_TASK_ID

echo "BEGIN DOWNLOAD"
datasets download genome accession --dehydrated --inputfile mydataset_${FILE}.txt --filename batch_${FILE}.zip

echo "BEGIN UNZIP"
unzip batch_${FILE}.zip -d batch${FILE}

echo "REHYDRATE"
datasets rehydrate --directory batch${FILE}

echo "CHECK MD5"
cd batch${FILE}
md5sum -c md5sum.txt > ../checks${FILE}.out

echo "END"
date