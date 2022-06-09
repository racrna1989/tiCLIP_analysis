#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

#wd

WD="/home/racrna/4_xiCLIP/"
cd $WD

##4_Map the snoRNA to the read

#then join the mature snoRNA file with the extended/mapped read. Then calculate the relative distance to the 3'end of the mature snoRNA
#if on negative take the left most coordinate from the snoRNA ($9) , from the right most of the read ($4)
#also notice it need a uniq

TMP1="$(mktemp racrna.XXXXXX.tmp)"

OUTDIR="10_snoRNA_read_coverage/3_map_snoRNA_to_read/"

INDIR="10_snoRNA_read_coverage/2_5end_and_3end_read_joined/" 

ANNODOWNSTREAMSNO="annotationFiles/snoRNAs.GRCh38andrefGene.maturetodowntreamexon.bed"
ANNOMATURESNO="annotationFiles/snoRNAs.GRCh38andrefGene.mature.bed"

mkdir -p $OUTDIR 

for EXTENDEDREAD in ${INDIR}*extended.bed;

do

ID=$(basename ${EXTENDEDREAD} | sed 's/\./\t/g' | cut -f 1)
#add step to exclude reads over 500nt
#step on is to map the snoRNA to the extended read. Delete the additional ":::5endTo3ss" and sort by snoRNA. Save as temp
bedtools intersect -s -loj \
-a <( sort -k1,1 -k2,2n ${EXTENDEDREAD} | awk '{OFS="\t"} $3-$2 < 501 {print $0}') \
-b <( sort -k1,1 -k2,2n ${ANNODOWNSTREAMSNO} ) | \
sed 's/:::5endTo3ss//g' | \
awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$10}' | sort -k7,7 > ${TMP1}


#if on + strand take the left most coordiates from the read ($3), from the right most of the snoRNA ($10)

join -1 7 -2 4 ${TMP1} <(sort -k4,4 ${ANNOMATURESNO} ) | \
sed 's/ /\t/g' | \
awk '{OFS="\t"}{
if($7 == "+"){
print $2,$3,$4,$5,$6,$7,($3-$10),($4-$10),$1} 
else if ($7=="-") 
{print  $2,$3,$4,$5,$6,$7,($9-$4),($9-$3),$1}
}'> ${OUTDIR}${ID}".relDistTo3endofmaturesnoRNA.extendedRead.bed"

rm ${TMP1}

done