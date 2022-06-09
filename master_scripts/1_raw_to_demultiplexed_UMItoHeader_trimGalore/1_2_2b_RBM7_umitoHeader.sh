#!/bin/bash
#SBATCH --mem-per-cpu 8G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP
###works!### 
#100419 - removes UMI from 9nt barcode at 5' end of read 1 (5nt) and 5' end of read 2 (1nt) and appends it to read 1 and 2 header.       

WD="/home/racrna/4_xiCLIP/"
cd $WD

inDir="2_demultiplexing/2_RT_demultiplexed_merged/"

for r1 in $inDir*NNN*RBM7*R1.fastq
do

r2=${r1/.R1.fastq/.R2.fastq}
umi_tools extract -I $r1 \
--bc-pattern=NNNXXXXNN \
--bc-pattern2=XXXXXXXXN \
--read2-in=$r2 \
--stdout=${r1/.fastq/.processed.fastq} \
--read2-out=${r2/.fastq/.processed.fastq} \
--log extract.log 

if [ -s "${r1/.fastq/.processed.fastq}" ]
then
rm $r1 $r2
else
echo "file exists but is empty" 
fi

done

