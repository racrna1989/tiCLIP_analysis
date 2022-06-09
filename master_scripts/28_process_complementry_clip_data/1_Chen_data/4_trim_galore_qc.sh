#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk


source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

WD="/home/racrna/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/1_Chen_ALYREF_iCLIP/"
cd $WD

INDIR="2_umiToHeader/"
OUTDIR="4_post_trim_galore/"

mkdir -p ${OUTDIR}

for INFQ in ${INDIR}*.UmiToHeader.fastq
do

trim_galore \
--stringency 4 \
--phred33 \
--fastqc \
--fastqc_args "--outdir 5_post_trim_galore_fastqc/" \
-a TGAGATCGGAAGAGCGG \
--clip_R1 4 \
--output_dir ${OUTDIR} \
${INFQ}

done


####unique RNA 3' adenylated adapter, and corrosponding iCLIP RT primer are different from that used in the standard iCLIP experiments. 
