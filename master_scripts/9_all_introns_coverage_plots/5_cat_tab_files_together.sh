#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/8_rna_binding_maps/1_intron_exon_junction_coverage/"
cd $WD

for INFILE in */*_biotype.tab; do 

OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 3 ) 


head -n 1 $INFILE > head.tmp

awk '{OFS="\t"} NR > 1 {print $0}' $INFILE >> "xiCLIP_all_"$OUTNAME".tmp"

done

cat head.tmp "xiCLIP_all_"$OUTNAME".tmp" > "xiCLIP_all_"$OUTNAME".counts"

rm head.tmp "xiCLIP_all_"$OUTNAME".tmp"


###

for INFILE in */*_sense.tab; do 

OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 3 ) 


head -n 1 $INFILE > head.tmp

awk '{OFS="\t"} NR > 1 {print $0}' $INFILE >> "xiCLIP_all_"$OUTNAME".tmp"

done

cat head.tmp "xiCLIP_all_"$OUTNAME".tmp" > "xiCLIP_all_"$OUTNAME".counts"

rm head.tmp "xiCLIP_all_"$OUTNAME".tmp"

