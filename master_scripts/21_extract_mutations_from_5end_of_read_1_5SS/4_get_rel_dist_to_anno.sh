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

ANNOBED="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.5SS.bed"
INDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_5end_of_read1_bed_overlaps/"
OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/6_rel_dist_from_annotation/" 
RANGE="10"
RANGE_NEG=$((RANGE *-1))
mkdir -p $OUTDIR 

for INBED in ${INDIR}*"5endofread1.bed";
do
ID=$(basename $INBED | cut -d "." -f 1 )

#bedtools intersect with loj to join bed files 

	
	echo "starting analysis for $ID"
	
	#intersect every read with closes annotation and print distance from read. antisense reads require fixing 
	
	bedtools closest \
	-a <( sort -k1,1 -k2,2n ${INBED} | sed 's/^chr//g' | cut -f 1-6 ) \
	-b <( sort -k1,1 -k2,2n ${ANNOBED} | sed 's/^chr//g' ) -D b | \
	
	#only print those reads within range 	
	awk -v RANGE=${RANGE} -v RANGE_NEG=${RANGE_NEG} '$NF >= RANGE_NEG && $NF <= RANGE' > ${OUTDIR}${ID}".rel_dist_to_anno.tab"

done
