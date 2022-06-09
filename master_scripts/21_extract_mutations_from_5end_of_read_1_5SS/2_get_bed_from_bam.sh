#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

##################
#
#Purpose - To extract read1 and ultimately extract the mutational profile of the reads 
#290620
#
##################


WD="/home/racrna/4_xiCLIP/"

cd $WD

INDIR="3_mapping/2_QC_mapped/"
OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/"

mkdir -p $OUTDIR

#filenames and dirs req


ANNOBED="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.5SS.bed"


OUTHEADERDIR=${OUTDIR}"1_samHeaders/"
OUTBAMDIR=${OUTDIR}"2_read1_bam_overlaps/"
OUTBEDDIR=${OUTDIR}"3_5end_of_read1_bed_overlaps/"


#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#make dirs if not already present
mkdir -p ${OUTHEADERDIR} ${OUTBAMDIR} ${OUTBEDDIR}


for INBAM in ${INDIR}*".SoUmiDedupRemSec.bam"; 
do 

#edit filename to save ID
ID=$(basename ${INBAM} | cut -d "." -f 1)

#1_save header from original read1 bamfile.
samtools view -f 64 -H ${INBAM} > ${OUTHEADERDIR}${ID}".SoUmiDedupRemSec.header" 


#mktmp file
TMPFILE_1=${OUTDIR}"tmp_1.tmp"


#2_extract read1 of reads that overlap bam region file. Bam file region is annotation extended by 10nt either side. Grep reads, cat header and extracted reads
samtools view -f 64 -L <( awk '{OFS="\t"}{$2=($2-10); $3=($3+10); print $0}' ${ANNOBED}) ${INBAM} | sort -k1,1 - > ${TMPFILE_1} 

cat ${OUTHEADERDIR}${ID}".SoUmiDedupRemSec.header" $TMPFILE_1 | samtools view -b | samtools sort - > ${OUTBAMDIR}${ID}".overlapping_annotation.read1.bam"

rm ${TMPFILE_1} 

#3_make bed file and convert it to cover only first nt of read
bedtools bamtobed -cigar -i ${OUTBAMDIR}${ID}".overlapping_annotation.read1.bam" | awk '{OFS="\t"}{
if($6 == "+"){print $1,$2,$2+1,$4,$5,$6,$7}
else if($6 == "-"){print $1,$3-1,$3,$4,$5,$6,$7}
}'> ${OUTBEDDIR}${ID}".overlapping_annotation.5endofread1.tmp.bed"

#get MD tag and paste next to bedfile

paste ${OUTBEDDIR}${ID}".overlapping_annotation.5endofread1.tmp.bed" <(samtools view ${OUTBAMDIR}${ID}".overlapping_annotation.read1.bam" | sed 's/.*MD:Z://g' | awk '{OFS="\t"}{print $1}') > ${OUTBEDDIR}${ID}".overlapping_annotation.5endofread1.bed"


rm ${OUTBEDDIR}${ID}".overlapping_annotation.5endofread1.tmp.bed"

done


rm ${OUTBEDDIR}"*val.overlapping*"

