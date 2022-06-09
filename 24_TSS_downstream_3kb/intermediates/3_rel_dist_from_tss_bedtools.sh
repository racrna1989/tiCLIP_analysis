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

INDIR="/home/racrna/4_xiCLIP/8_rna_binding_maps/6_PROMPTS/5primepos/1_overlaps/"
OUTDIR="/home/racrna/4_xiCLIP/8_rna_binding_maps/6_PROMPTS/5primepos/2_rel_dist_to_TSS/"



#constructing output suffix for filename. 
OUTFILENAME=$(basename $ANNO | sed 's/.bed//g')

#mk output dir

mkdir -p ${OUTDIR}

#preform analysis

for INBED in ${INDIR}*.bed;
	do
	#make variables
	
	ID=$(basename ${INBED} | cut -d "." -f 1 )

	echo "starting analysis for $ID"
	
	bedtools closest \
	-a <( sort -k1,1 -k2,2n ${INBED} | sed 's/^chr//g' ) \
	-b <( sort -k1,1 -k2,2n ${ANNO} | sed 's/^chr//g' ) -D b | \
	awk -v ID=${ID} '{OFS="\t"}{
	if($6 != $12){
	$NF=($NF *= -1); print ID"_antisense",$0}
	else if ($6 == $12){
	print ID"_sense",$0}}' | awk '$NF >= -3002 && $NF <= 3002' > ${OUTDIR}${ID}"."${OUTFILENAME}".dist_to_TSS.tab"
	
	
done




