#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/8_rna_binding_maps/2_intron_exon_junction_coverage_retained_introns/"
cd $WD

#5end

for INFILE in */*no_gene_norm.tab; do 

OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 2 ) 

head -n 1 $INFILE > head.tmp

awk '{OFS="\t"} NR > 1 {print $0}' $INFILE >> "xiCLIP_all_"$OUTNAME".tmp"

done

cat head.tmp "xiCLIP_all_"$OUTNAME".tmp" > "xiCLIP_all_"$OUTNAME".tab"

rm head.tmp "xiCLIP_all_"$OUTNAME".tmp"

