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


#variables that contain filepaths for input fq, output bam and output qc'ed bam
FQINDIR="1_fastqgz/"
OUTDIR="3_fastqgz_all/"

mkdir -p ${OUTDIR}


#cat all sequencing runs together for each 

PROTEIN_LIST='Alyref Chtop Cntrl Nxf1'
REP_LIST='R1 R2 R3'


INDIR="1_fastqgz/"

for P in ${PROTEIN_LIST}; do
	for R in ${REP_LIST}; do
	
	#cat all files from same replicate 
	zcat ${INDIR}/*${P}*${R}* | sed 's/ /_/g' > ${OUTDIR}${P}"_FLAG_"${R}"_all.fastq"
	gzip ${OUTDIR}${P}"_FLAG_"${R}"_all.fastq"

	done
done

