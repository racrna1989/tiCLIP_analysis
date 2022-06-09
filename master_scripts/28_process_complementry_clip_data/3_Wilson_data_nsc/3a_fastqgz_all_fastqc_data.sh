#!/bin/bash
#SBATCH --mem-per-cpu 8G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

WD="/home/racrna/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/3_Wilson_EJC_iCLIP_softCLIP/"

cd ${WD}

INDIR="3_fastqgz_all/"
OUTDIR="3a_fastqz_all_fastqc/"

mkdir -p ${OUTDIR}

# preform fastqc 


for INFQ in ${INDIR}*.gz;
do

fastqc ${INFQ} \
--outdir ${OUTDIR} \
-t 8 \
-q 

done
