#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

wd="/home/racrna/4_xiCLIP/"
cd $wd


#making updated snRNA annotation file. 

#step make annotation file for 10 nt up and down from 5SS and BP

BP_ANNO="annotationFiles/BP-annotated-metadata.hg38.bed"
INTRON_ANNO="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.bed"

#get5primeFlank.bash $BP_ANNO 10 10 | sed 's/^chr//g' - > ${BP_ANNO/.bed/.5.bed}
get5primeFlank.bash $INTRON_ANNO 0 0 | sed 's/^chr//g' - > ${INTRON_ANNO/.bed/.5SS.bed}



