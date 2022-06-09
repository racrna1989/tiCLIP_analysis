#!/bin/bash
#SBATCH --mem-per-cpu 40G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP
###works!### 
#100419 - removes UMI from 9nt barcode at 5' end of read 1 (5nt) and 5' end of read 2 (1nt) and appends it to read 1 and 2 header.       



WD="/home/racrna/xiCLIP/faststorage/Projects/xiCLIP_data/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/1_Chen_ALYREF_iCLIP"

cd ${WD}

INDIR="1_fastqgz/"
OUTDIR="2_umiToHeader/"

mkdir -p ${OUTDIR}

for FQIN in ${INDIR}*fastq.gz; 
do 

OUTID=$(basename $FQIN | sed 's/.gz//g' )
FQOUT=$(echo ${OUTDIR}${OUTID})

zcat ${FQIN} | paste - - - - | cut -d " " -f 2 | sed 's/\t/\n/g' | \
umi_tools extract \
--bc-pattern=NNNXXXXNN \
--stdout=${FQOUT/.fastq/.UmiToHeader.fastq} \
--log=${FQIN/.fastq/.UmiToHeader.log}


done
