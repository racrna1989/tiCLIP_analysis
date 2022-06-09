#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk



##################
#
#Purpose is to split bam files in to splice ($6 contains an N). and not split #Purpose - To extract read1 and ultimately extract the mutational profile of the reads 
#171220
#
##################

#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#set wd
WD="/home/racrna/4_xiCLIP/"
cd ${WD}

#filenames and dirs req


INDIR="3_mapping/2_QC_mapped/"
MASTEROUTDIR="12_extract_splice_and_non_spliced_reads_close_to_exon_junctions/"
OUTDIR=${MASTEROUTDIR}"1_spliced_and_not_spliced_bams/"

mkdir -p ${MASTEROUTDIR} ${OUTDIR}


for INBAM in ${INDIR}*SoUmiDedupRemSec.bam
do
ID=$(basename ${INBAM} | cut -d "." -f 1 )

#make spliced and not spliced files for analysis for read 1

#step 1 is to extract read1 of properly paired mapped reads that are split or not split (N in cigar string)  
samtools view -@ 8 -f 67 $INBAM | awk 'IFS=OFS="\t" $6~"N" {print $0}' > ${OUTDIR}${ID}.spliced.sam
samtools view -@ 8 -f 67 $INBAM | awk 'IFS=OFS="\t" $6!~"N" {print $0}' > ${OUTDIR}${ID}.notspliced.sam

#step 2 add header back to sam and convert back to sorted bam.
samtools view -H $INBAM | cat - ${ID}.spliced.sam | samtools sort - | samtools view -b - > ${OUTDIR}${ID}.spliced.read1.bam
samtools view -H $INBAM | cat - ${ID}.notspliced.sam | samtools sort - | samtools view -b - > ${OUTDIR}${ID}.notspliced.read1.bam
rm *.sam
wait

#make spliced and not spliced files for analysis for read 2

#step 1 is to extract read2 of properly paired mapped reads that are split or not split (N in cigar string)  
samtools view -@ 8 -f 131 $INBAM | awk 'IFS=OFS="\t" $6~"N" {print $0}' > ${OUTDIR}${ID}.spliced.sam
samtools view -@ 8 -f 131 $INBAM | awk 'IFS=OFS="\t" $6!~"N" {print $0}' > ${OUTDIR}${ID}.notspliced.sam

#step 2 add header back to sam and convert back to sorted bam.

samtools view -H $INBAM | cat - ${ID}.spliced.sam | samtools sort - | samtools view -b - > ${OUTDIR}${ID}.spliced.read2.bam
samtools view -H $INBAM | cat - ${ID}.notspliced.sam | samtools sort - | samtools view -b - > ${OUTDIR}${ID}.notspliced.read2.bam
rm *.sam
wait

done
