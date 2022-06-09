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

INDIR="/home/racrna/4_xiCLIP/8_rna_binding_maps/6_PROMPTS/3endOfRead2/1_overlaps/"
OUTDIR="/home/racrna/4_xiCLIP/8_rna_binding_maps/6_PROMPTS/3endOfRead2/2_rel_dist_to_TSS/"



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
	
	bedtools intersect -loj \
	-a <( sort -k1,1 -k2,2n ${INBED} | sed 's/^chr//g' ) \
	-b <( awk '{OFS="\t"}{print $1,$2-3001,$3+3001,$4,$5,$6}' ${ANNO} | sort -k1,1 -k2,2n | sed 's/^chr//g' ) | sort -k4,4 -k10,10 | \
	awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$10}' > ${OUTDIR}${ID}"."${OUTFILENAME}"_joined_with_PROMPT.counts"
	
	
	
	#join files based on common ID.
	sort -k7,7 ${OUTDIR}${ID}"."${OUTFILENAME}"_joined_with_PROMPT.counts"| join -1 7 -2 4 - <(sort -k4,4 ${ANNO} ) | tr " " "\t" | sort -k1,1 -k2,2n | \
	#calculation based on strand to determine distance of 5' pos from annotation
	awk '{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1,$11,$12}' | \
	awk -v ID=${ID} '{OFS="\t"}{
	if($6 == "+" && $6 == $12){
	print ID"_sense",$1,$2,$3,$4,$5,$6,$10,$11,$12,($2-$8)} 
	else if($6 == "+" && $6 != $12){
	print ID"_antisense",$1,$2,$3,$4,$5,$6,$10,$11,$12,($2-$8)} 
	else if ($6=="-" && $6 == $12){
	print ID"_sense",$1,$2,$3,$4,$5,$6,$10,$11,$12,($8-$2)}
	else if ($6=="-" && $6 != $12){
	print ID"_antisense",$1,$2,$3,$4,$5,$6,$10,$11,$12,($8-$2)}
	}' | sort -k5,5 -k8,8 | uniq > ${OUTDIR}${ID}"."${OUTFILENAME}".counts"
	#
	#echo "completed coverage counts for $ID"
	#
done




