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

SCRIPT_TSS="master_scripts/24_TSS_downstream_3kb/rel_dist_TSS.sh"
SCRIPT_TTS="master_scripts/24_TSS_downstream_3kb/rel_dist_TTS.sh"
OUTDIR_MASTER="/home/racrna/4_xiCLIP/8_rna_binding_maps/7_all_genes_TTS_TSS/"
RANGE="3002"
ANNO="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed"


mkdir -p ${OUTDIR_MASTER}

for READTYPE in 5primepos 3endOfRead2 ; 

	do
	
	sbatch ${SCRIPT_TSS} ${OUTDIR_MASTER} ${READTYPE} ${RANGE} ${ANNO}
	sbatch ${SCRIPT_TTS} ${OUTDIR_MASTER} ${READTYPE} ${RANGE} ${ANNO}

done
