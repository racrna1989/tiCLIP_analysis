#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 10G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk


WD="/home/racrna/4_xiCLIP/"

cd $WD

INDIR_MISMATCHES="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_getMutProfile/4_read_mismatches/"
INDIR_RELDIST="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/2_dist_from_annotation/"
ANNOBED="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.5SS.bed"

OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/4_relDist_to_anno_and_mismatches/"

mkdir -p $OUTDIR


# join rel dist and mutation file based on read id	
for INTAB in ${INDIR_MISMATCHES}*_wrangled.tab; 
do

ID=$(basename $INTAB | cut -d "." -f 1)

TMPFILE_1=${OUTDIR}${ID}".relDist.tmp"
TMPFILE_2=${OUTDIR}${ID}".mismatch.tmp"


sed "s|/1\t|\t|g" ${INDIR_RELDIST}${ID}".relDistOf5primeposToAnno.bed" | cut -f 1-9 | sort -k4,4 > ${TMPFILE_1}

sort -k1,1 ${INTAB} > ${TMPFILE_2}

join -1 4 -2 1 ${TMPFILE_1} ${TMPFILE_2} | tr " " "\t" > ${OUTDIR}${ID}".mutations_and_rel_dist_from_annotation.tab"

rm ${TMPFILE_1} ${TMPFILE_2}

done
#remove final file if it exits 
rm -f ${OUTDIR}"xiCLIP_all.mutations_and_rel_dist_from_annotation.tab"

#cat all files together and append name to first feild

for INTAB in ${OUTDIR}*".mutations_and_rel_dist_from_annotation.tab"
do
ID=$(basename $INTAB | cut -d "." -f 1)

awk -v ID=$ID '{OFS="\t"}{print ID,$0}' $INTAB >> ${OUTDIR}"xiCLIP_all.mutations_and_rel_dist_from_annotation.tab"

done