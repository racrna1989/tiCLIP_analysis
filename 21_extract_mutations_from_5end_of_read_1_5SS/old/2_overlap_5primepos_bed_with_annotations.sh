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
SPLICESITE5SS_ANNO="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.5SS.bed"
#BP_ANNO="annotationFiles/BP-annotated-metadata.hg38.updown10nt.bed"

#OUTDIR_BP="/home/racrna/4_xiCLIP/11_extract_mutations_from_5end_of_read_1/1_overlap_5primepos_reads_with_annotation/branchpoints/"
OUTDIR_5SS="/home/racrna/4_xiCLIP/11_extract_mutations_from_5end_of_read_1/1_5splice_sites/1_overlap_5primepos_reads_with_annotation/"

mkdir -p $OUTDIR_BP $OUTDIR_5SS

INDIR="5_bedGraphs_and_derivatives/5primepos/"

for INBED in ${INDIR}*.bed 
do
NAME=$(basename $INBED | sed 's/\./\t/g' | cut -f 1 )

# add 10nt window up and down
bedtools intersect -s -a ${INBED} -b <( awk '{OFS="\t"}{$2=($2-10); $3=($3+10); print $0}' ${SPLICESITE5SS_ANNO} ) > ${OUTDIR_5SS}${NAME}".5SS_overlap.bed" &

wait

done


