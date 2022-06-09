#!/bin/bash

# Master variables required for script 

OUTDIR_MASTER=$1
READTYPE=$2
RANGE=$3
ANNO=$4

#1 make TSS annotation


TTSANNO=${ANNO/.bed/.TTS.bed}
if [ -s "$file" ]

then 
   echo " file exists and is not empty "
else
   echo " file does not exist, or is empty "
      
	awk '{OFS="\t"}{
	if ($6 == "+"){
	$3=$2+1; print $0}
	else if ($6 == "-"){
	$2=$3-1; print $0}}' ${ANNO} > ${TTSANNO}
	   
fi


#2 overlap 5'pos with extended snoRNA annotation 

INDIR="5_bedGraphs_and_derivatives/5primepos/"
OUTDIR=${OUTDIR_MASTER}${READTYPE}"/1_overlaps/"
OUTSUFFIX=$(basename ${TTSANNO} | sed 's/.bed//g')

mkdir -p $OUTDIR

for INBED in ${INDIR}*.bed 
	do
	NAME=$(basename $INBED | cut -d "." -f 1 )

	# add 3000nt window up and down. also capture reads on opposite strand
	bedtools intersect \
	-a ${INBED} \
	-b <( awk '{OFS="\t"}{print $1,$2-3001,$2+3001,$4,$5,$6}' ${TTSANNO} ) \
	> ${OUTDIR}${NAME}"."${OUTSUFFIX}".bed"

done


#3 calculate rel distance

INDIR=${OUTDIR}
OUTDIR=${OUTDIR_MASTER}${READTYPE}"/2_rel_dist_to_TTS/"
OUTFILENAME=$(basename ${TTSANNO} | sed 's/TTS.bed/dist_to_TTS.tab/g')
	
RANGE_NEG=$((RANGE *-1))
	
mkdir -p ${OUTDIR}

for INBED in ${INDIR}*.bed;
	do

	
	ID=$(basename ${INBED} | cut -d "." -f 1 )

	echo "starting analysis for $ID"
	
	#intersect every read with closes annotation and print distance from read. antisense reads require fixing 
	
	bedtools closest \
	-a <( sort -k1,1 -k2,2n ${INBED} | sed 's/^chr//g' ) \
	-b <( sort -k1,1 -k2,2n ${TTSANNO} | sed 's/^chr//g' ) -D b | \
	
	#fix distance upstream downstream for antisense
	
	awk -v ID=${ID} -v READTYPE=${READTYPE} '{OFS="\t"}{
	if($6 != $12){
	print ID"_antisense_TTS",$0}
	else if ($6 == $12){
	print ID"_sense_TTS",$0}
	}' | \
	
	#only print those reads within range 	
	awk -v RANGE_NEG=$RANGE_NEG -v RANGE=${RANGE} '$NF >= RANGE_NEG && $NF <= RANGE' > ${OUTDIR}${ID}"."${OUTFILENAME}

done


#4 cat together files

rm ${OUTDIR}"xiCLIP_all_"${READTYPE}"."${OUTFILENAME}

for INTAB in ${OUTDIR}*.tab ;
	do
	
	awk '{OFS="\t"}{print $1,$8,$9,$10,$11,$12,$13,$14}' ${INTAB} >> ${OUTDIR}"xiCLIP_all_"${READTYPE}"."${OUTFILENAME}

done