#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

WD="/home/racrna/4_xiCLIP/"
cd $WD


inBamDir="3_mapping/2_QC_mapped/"
inBed="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.77.cleaned.sorted.2.gene.bed"

outDir="4_geneCount/"
outCountDirAllGenes="4_geneCount/allGenescount/"
outrRNACountDir="4_geneCount/rRNA/"

mkdir -p $outrRNACountDir $outCountDirAllGenes

for inBam in $inBamDir*CBP80*"SoUmiDedupRemSec.bam"
do
name=$(basename $inBam | sed 's/.SoUmiDedupRemSec.bam//g')

## count overlap with all genes

bedtools intersect -c -s -a <(sort -k1,1 -k2,2n $inBed) -b $inBam > $outCountDirAllGenes$name".allGenescount.tab" 

##count overlap with rRNA genes

bedtools intersect -c -s -a <(grep rRNA $inBed | sort -k1,1 -k2,2n) -b $inBam > $outrRNACountDir$name".rRNAcount.tab" 

## count total number of reads (these are properly paired mapped reads)

count=$(samtools view -c $inBam) 

echo -e $name"\t"$count >> 4_geneCount/allCounts.txt

done 
