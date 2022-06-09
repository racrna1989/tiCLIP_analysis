#!/bin/bash
#SBATCH --mem-per-cpu 8G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH -t 2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/"
cd $WD

subDir="2_demultiplexing/"
inDir="2_demultiplexing/1_RT_demultiplexed/"
outDir="2_demultiplexing/2_RT_demultiplexed_merged/"

mkdir -p $subDir $inDir $outDir

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyCRAC-tools

#generate list of samples, by removing non-redunant characters from identifier#

ls $inDir*.fastq | sed 's|/|\t|g' | cut -f 3 | cut -c 15- | rev | cut -c 7- | rev | sort | uniq > $subDir"sample.list"

for f in RBM7 ALYREF CBP20 CBP80 PHAX 
do
cat $subDir"sample.list" | sed 's/_/\t/g' | cut -f 2,3 | sed 's/\t/_/g' | grep $f > $subDir$f"_sample.list"
done

#read through the samples lis and use identifer in ls function, plus the read number, and cat these files. file saved with extention of merged.R?.fastq.gz

for SAMPLENAME in $(cat $subDir"sample.list"); do ls $inDir*R1*$SAMPLENAME*| xargs cat > $outDir$SAMPLENAME".merged.R1.fastq" ; done &
for SAMPLENAME in $(cat $subDir"sample.list"); do ls $inDir*R2*$SAMPLENAME*| xargs cat > $outDir$SAMPLENAME".merged.R2.fastq" ; done &

wait
