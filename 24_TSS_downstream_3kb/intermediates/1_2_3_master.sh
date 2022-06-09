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

# Master variables required for script 

OUTDIR_MASTER="/home/racrna/4_xiCLIP/8_rna_binding_maps/6_PROMPTS/"
READTYPE="5primepos"
RANGE="3002"
ANNO="annotationFiles/PROMPTs_and_HeLa_PROMPTs.bed"

#1 make TSS annotation

TSSANNO=${ANNO/.bed/.TSS.bed}

awk '{OFS="\t"}{
if ($6 == "+"){
$3=$2+1; print $0}
else if ($6 == "-"){
$2=$3-1; print $0}}' ${ANNO} > ${TSSANNO}

#2 overlap 5'pos with extended snoRNA annotation 

INDIR="5_bedGraphs_and_derivatives/5primepos/"
OUTDIR=${OUTDIR_MASTER}${READTYPE}"/1_overlaps/"
OUTSUFFIX=$(basename ${TSSANNO} | sed 's/.bed//g')

mkdir -p $OUTDIR

for INBED in ${INDIR}*.bed 
	do
	NAME=$(basename $INBED | cut -d "." -f 1 )

	# add 3000nt window up and down. also capture reads on opposite strand
	bedtools intersect \
	-a ${INBED} \
	-b <( awk '{OFS="\t"}{print $1,$2-3001,$2+3001,$4,$5,$6}' ${TSSANNO} ) \
	> ${OUTDIR}${NAME}"."${OUTSUFFIX}".bed"

done


#3 calculate rel distance

INDIR=${OUTDIR}
OUTDIR=${OUTDIR_MASTER}${READTYPE}"/2_rel_dist_to_TSS/"
OUTFILENAME=$(basename ${TSSANNO} | sed 's/TSS.bed/dist_to_TSS.tab/g')
mkdir -p ${OUTDIR}

for INBED in ${INDIR}*.bed;
	do

	
	ID=$(basename ${INBED} | cut -d "." -f 1 )

	echo "starting analysis for $ID"
	
	#intersect every read with closes annotation and print distance from read. antisense reads require fixing 
	
	bedtools closest \
	-a <( sort -k1,1 -k2,2n ${INBED} | sed 's/^chr//g' ) \
	-b <( sort -k1,1 -k2,2n ${ANNO} | sed 's/^chr//g' ) -D b | \
	
	#fix distance upstream downstream for antisense
	
	awk -v ID=${ID} '{OFS="\t"}{
	if($6 != $12){
	$NF=($NF *= -1); print ID"_antisense",$0}
	else if ($6 == $12){
	print ID"_sense",$0}}' | \
	
	#only print those reads within range 
	
	RANGE_NEG=$((RANGE *-1))
	
	awk -v RANGE_NEG=$RANGE_NEG -v RANGE=${RANGE} '$NF >= RANGE_NEG && $NF <= RANGE ' > ${OUTDIR}${ID}"."${OUTFILENAME}
	
	
done

