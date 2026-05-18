#!/usr/bin/env bash
#SBATCH -J tax_grab
#SBATCH -c 1
#SBATCH --output=%x_%j.out

echo "START"
date

cd $MYPATH/entropy/
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

echo "END"
date
