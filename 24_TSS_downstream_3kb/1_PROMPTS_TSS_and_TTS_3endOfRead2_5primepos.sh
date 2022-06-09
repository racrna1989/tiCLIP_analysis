#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk
#set working directory

wd="/home/racrna/4_xiCLIP/"
cd $wd

SCRIPT="master_scripts/24_TSS_downstream_3kb/rel_dist_TSS.sh"
OUTDIR_MASTER="/home/racrna/4_xiCLIP/8_rna_binding_maps/6_PROMPTS_TSS/"
RANGE="3002"
ANNO="annotationFiles/PROMPTs_and_HeLa_PROMPTs.bed"

for READTYPE in 5primepos 3endOfRead2 ; 

	do
	
	sbatch ${SCRIPT} ${OUTDIR_MASTER} ${READTYPE} ${RANGE} ${ANNO}

done