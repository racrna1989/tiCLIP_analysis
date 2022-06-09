#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

#wd

wd="/home/racrna/4_xiCLIP/"
cd $wd

#overlap 5'pos with extended snoRNA annotation 
INANNO="annotationFiles/snoRNAs.GRCh38andrefGene.maturetodowntreamexon.bed"
OUTDIR="/home/racrna/4_xiCLIP/10_snoRNA_read_coverage/1_intersect_with_snoRNA_and_3prime_extention/"
INDIR="5_bedGraphs_and_derivatives/5primepos/"

mkdir -p $OUTDIR

for INBED in ${INDIR}*.bed 
do
NAME=$(basename $INBED | sed 's/\./\t/g' | cut -f 1 )

bedtools intersect -s -a ${INBED} -b ${INANNO} > ${OUTDIR}${NAME}.overlap.bed

done

