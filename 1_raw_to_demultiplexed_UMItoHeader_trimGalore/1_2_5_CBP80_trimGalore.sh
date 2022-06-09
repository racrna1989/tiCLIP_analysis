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

WD="/home/racrna/4_xiCLIP/"
cd $WD

inDir="2_demultiplexing/3_RT_L3_demultiplexed/"
outDir="2_demultiplexing/4_RT_L3_demultiplexed_trimGalore/"
mkdir -p $outDir

for R1 in $inDir*CBP80*R1.processed*.fastq
do

R2=$(echo $R1 | sed 's/R1\./R2\./g;s/_1\.fastq/_2\.fastq/g')

trim_galore \
--stringency 4 \
--phred33 \
--fastqc \
--illumina \
-a2 AGATCGGAAGAGC \
--clip_R1 4 \
--clip_R2 8 \
--three_prime_clip_R1 9 \
--three_prime_clip_R2 9 \
--trim1 \
--output_dir $outDir \
--paired $R1 $R2 

done



