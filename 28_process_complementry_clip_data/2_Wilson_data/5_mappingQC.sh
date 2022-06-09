#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 10
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#make new working directory and cd to it
WD="/home/racrna/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/2_Wilson_EJC_iCLIP/"
cd $WD

#variables that contain filepaths for input fq, output bam and output qc'ed bam

INBAMDIR="4_mapped/"
OUTMAPQCDIR="5_QC_mapped/"

mkdir -p ${OUTMAPQCDIR}

#mapping with no_softClipping
REFERENCE=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa
GENOME=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned
GTF=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.77.cleaned.sorted.2.gtf
SS=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.77.cleaned.sorted.2.ss

for INBAM in ${INBAMDIR}*".sam";
do

ID=$(basename $INBAM | sed 's/.sam//g')

##sort sam to bam
samtools import ${REFERENCE} ${OUTBAMDIR}${ID}".sam" - | samtools sort -@ 10 - -o ${OUTMAPQCDIR}${ID}".So.bam" 

##index bam for dedup
samtools index -@ 10 ${OUTMAPQCDIR}${ID}".So.bam" 

##remove duplicates based on umi in header 
umi_tools dedup -I ${OUTMAPQCDIR}${ID}".So.bam" --output-stats=${OUTMAPQCDIR}${ID}".dedupStats." -S ${OUTMAPQCDIR}${ID}".SoUmiDedup.bam" ;

##remove secondary alignments and nonmapped reads
samtools view -b -F 260 ${OUTMAPQCDIR}${ID}".SoUmiDedup.bam" > ${OUTMAPQCDIR}${ID}".SoUmiDedupRemSec.bam"

done
