#!/bin/bash
#SBATCH --mem-per-cpu 8G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH -t 2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

#preform second demultiplex step

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyCRAC-tools

WD="/home/racrna/4_xiCLIP/"
subDir="2_demultiplexing/"
inDir="2_demultiplexing/2_RT_demultiplexed_merged/"
outDir="2_demultiplexing/3_RT_L3_demultiplexed/"
bcDir="2_demultiplexing/barcodes/L3BC/"

cd $WD
mkdir -p $outDir 

for SAMPLENAME in $(cat $subDir"CBP20_sample.list"); do \
r1=$(ls $inDir*$SAMPLENAME*merged*R1*); \
r2=$(ls $inDir*$SAMPLENAME*merged*R2*); \
bc=$(ls $bcDir*$SAMPLENAME*); \

echo $r1
echo $r2
echo $bc


pyBarcodeFilter.py \
-f $r1 \
-r $r2 \
--file_type=fastq \
-b $bc \
-m 1 \
-k \
--search reverse 

wait 

echo "$SAMPLENAME complete" 

mv *$SAMPLENAME*.fastq $outDir
mv *$SAMPLENAME*.txt $outDir 

done

wait 


