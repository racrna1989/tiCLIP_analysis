#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/8_rna_binding_maps/4_intron_exon_junction_coverage_by_exon_context/"
cd $WD

for INFILE in */*.tab; do 

OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 2 ) 


head -n 1 $INFILE > "head."$OUTNAME".tmp"

awk '{OFS="\t"} NR > 1 {print $0}' $INFILE >> "xiCLIP_all."$OUTNAME".tmp"

done


for INFILE in *.tmp  ; do 

OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 2 ) 

cat "head."$OUTNAME".tmp" "xiCLIP_all."$OUTNAME".tmp" > "xiCLIP_all."$OUTNAME".tab"

#rm "head."$OUTNAME".tmp" "xiCLIP_all."$OUTNAME".tmp"

done

