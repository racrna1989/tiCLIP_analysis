#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk
#set working directory

wd="/home/racrna/4_xiCLIP/"
cd $wd

SCRIPT="master_scripts/25_miRNA-5p/get_rel_dist_from_anno.sh"
OUTDIR_MASTER="/home/racrna/4_xiCLIP/8_rna_binding_maps/8_miRNA-5p_5enda/"
RANGE="100"
ANNO="annotationFiles/miRNA_5p_5end.bed"


mkdir -p ${OUTDIR_MASTER}

for READTYPE in 5primepos 3endOfRead2 ; 

	do
	
	sbatch ${SCRIPT} ${OUTDIR_MASTER} ${READTYPE} ${RANGE} ${ANNO}

done
