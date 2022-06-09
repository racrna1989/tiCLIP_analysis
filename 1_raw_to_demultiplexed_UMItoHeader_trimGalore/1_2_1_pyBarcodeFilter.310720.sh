#!/bin/bash
#SBATCH --mem-per-cpu 8G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH --mail-user ross.cordiner@mbg.au.dk
#SBATCH --mail-type ALL
#SBATCH -t 24:00:00

WD="/home/racrna/4_xiCLIP/"
cd $WD

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyCRAC-tools

inDir="1_data/"
outDir="2_demultiplexing/1_RT_demultiplexed/"
bc="/home/racrna/4_xiCLIP/2_demultiplexing/barcodes/rtBC/rtBarcodes_xiCLIP-200704.tab"
mkdir -p $outDir


#read L001
pyBarcodeFilter.py \
-f 1_data/L001_R1.fastq.gz \
-r 1_data/L001_R2.fastq.gz \
--file_type=fastq.gz \
-b $bc \
-m 1 \
-k  &

#read L002
pyBarcodeFilter.py \
-f 1_data/L002_R1.fastq.gz \
-r 1_data/L002_R2.fastq.gz \
--file_type=fastq.gz \
-b $bc \
-m 1 \
-k &

#read L003
pyBarcodeFilter.py \
-f 1_data/L003_R1.fastq.gz \
-r 1_data/L003_R2.fastq.gz \
--file_type=fastq.gz \
-b $bc \
-m 1 \
-k &

#read L004
pyBarcodeFilter.py \
-f 1_data/L004_R1.fastq.gz \
-r 1_data/L004_R2.fastq.gz \
--file_type=fastq.gz \
-b $bc \
-m 1 \
-k &

wait

mv *.fastq $outDir
mv *.txt $outDir 




