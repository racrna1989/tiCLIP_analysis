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
ANNO="annotationFiles/PROMPTs_and_HeLa_PROMPTs.TSS.bed"

OUTDIR="/home/racrna/4_xiCLIP/8_rna_binding_maps/6_PROMPTS/5primepos/1_overlaps/"

INDIR="5_bedGraphs_and_derivatives/5primepos/"

mkdir -p $OUTDIR



for INBED in ${INDIR}*.bed 
do
NAME=$(basename $INBED | cut -d "." -f 1 )
OUTSUFFIX=$(basename $ANNO | sed 's/.bed//g')

# add 3000nt window up and down. also capture reads on opposite strand
bedtools intersect -a ${INBED} -b <( awk '{OFS="\t"}{print $1,$2-3001,$2+3001,$4,$5,$6}' ${ANNO} ) > ${OUTDIR}${NAME}"."${OUTSUFFIX}".bed"

done

