#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/8_rna_binding_maps/5_coverage_around_branch_points/"
cd $WD

for INFILE in */*.counts; do 

READNAME=$(dirname $INFILE | sed 's/-/\t/g' | cut -f 1 )
OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 4-7 | sed 's/\t/_/g') 

awk '{OFS="\t"}{print $0}' $INFILE >> "xiCLIP."$READNAME"."$OUTNAME".count"

done



