#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 10
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#make new working directory and cd to it
WD="/home/racrna/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/3_Wilson_EJC_iCLIP_softCLIP/"
cd $WD

INDIR="3_fastqgz_all/"
OUTBAMDIR="4_mapped/"

mkdir -p ${OUTBAMDIR}

#filepaths for mapping with no_softClipping
REFERENCE=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa
GENOME=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned
GTF=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.77.cleaned.sorted.2.gtf
SS=/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.77.cleaned.sorted.2.ss

for INFQGZ in ${INDIR}*".fastq.gz"; 
	do
	
	ID=$(basename ${INFQGZ} | sed 's/.fastq.gz//g')
	
	#hisat2 new parameters, no-softclip. 
	hisat2 \
	-t \
	--rna-strandness FR \
	--fr \
	--phred33 \
	-p 10 \
	--new-summary \
	--no-softclip \
	--known-splicesite-infile ${SS} \
	-x ${GENOME} \
	-U ${INFQGZ} \
	-S ${OUTBAMDIR}${ID}".sam"
	
done


