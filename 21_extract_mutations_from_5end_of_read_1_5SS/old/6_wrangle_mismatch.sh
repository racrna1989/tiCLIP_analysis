#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk


WD="/home/racrna/4_xiCLIP/"

cd $WD

OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_getMutProfile/4_read_mismatches/"
INDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_getMutProfile/4_read_mismatches/"
mkdir -p $OUTDIR



for INTAB in ${INDIR}*.read1_mismatch.tab ; 
do

awk '{OFS="\t"}{print $1,$2,$5,$6,$9,$10,$(NF-4)}' $INTAB > ${INTAB/.tab/_wrangled.tab}

done
