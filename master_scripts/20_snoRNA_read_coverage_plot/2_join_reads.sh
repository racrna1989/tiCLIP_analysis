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

INDIR="10_snoRNA_read_coverage/1_intersect_with_snoRNA_and_3prime_extention/"
OUTDIR="10_snoRNA_read_coverage/2_5end_and_3end_read_joined/" 
mkdir -p $OUTDIR 

INDIR3ENDOFREAD2="5_bedGraphs_and_derivatives/3endOfRead2/"



for INBED in ${INDIR}*overlap.bed
do
	#ID 
	ID=$(basename ${INBED} | sed 's/\./\t/g' | cut -f 1 )
	echo ${ID}
	#3endflank
	join -1 4 -2 4 \
	<(sed 's|/1||g' ${INBED} | sort -k4,4) \
	<(sed 's|/2||g' ${INDIR3ENDOFREAD2}${ID}*".bed" | sort -k4,4) | sed 's/ /\t/g' | \
	awk '{OFS="\t"}{
	if ($6 == "-" ){
	print $2, $8, $4, $1, $5, $6}
	else if ($6 == "+") {
	print $2, $3, $9, $1, $5, $6}
	}' > ${OUTDIR}${ID}".extended.bed"
	
done

