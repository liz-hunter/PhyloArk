#!/usr/bin/env bash
#SBATCH -J inspect
#SBATCH -n 5

module load kraken2/2.1.5

echo "START"
date

kraken2-inspect --db refseq_db > inspect/refseq_inspect.txt

echo "END"
date
