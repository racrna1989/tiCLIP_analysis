#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/8_rna_binding_maps/3_snoRNA_and_host_intron_coverage_at_ends_of_annotated_genes/"
cd $WD

#gene norm

#for INFILE in */*.tab; do 

#OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 2 ) 

#head -n 1 $INFILE > head.tmp

#awk '{OFS="\t"} NR > 1 {print $0}' $INFILE >> "xiCLIP_all_"$OUTNAME".tmp"

#done

#cat head.tmp "xiCLIP_all_"$OUTNAME".tmp" > "xiCLIP_all_"$OUTNAME".tab"

#rm head.tmp "xiCLIP_all_"$OUTNAME".tmp"

#no gene norm

for INFILE in */*no_gene_norm.tab; do 

OUTNAME=$(basename $INFILE | sed 's/\./\t/g' | cut -f 2 ) 

head -n 1 $INFILE > head.tmp

awk '{OFS="\t"} NR > 1 {print $0}' $INFILE >> "xiCLIP_all_"$OUTNAME".tmp"

done

cat head.tmp "xiCLIP_all_"$OUTNAME".tmp" > "xiCLIP_all_"$OUTNAME"_no_gene_norm.tab"

rm head.tmp "xiCLIP_all_"$OUTNAME".tmp"