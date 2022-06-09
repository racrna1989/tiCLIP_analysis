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


inBamDir="3_mapping/3_QC_mapped/"
inBed="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.77.cleaned.sorted.2.gene.bed"

outDir="4_geneCount/"
outCountDirAllGenes="4_geneCount/allGenescount/"
outrRNACountDir="4_geneCount/rRNA/"


#aggregate rRNA counts into one file.

for inTab in $outrRNACountDir*.tab ;

do
name=$(basename $inTab | sed 's/.rRNAcount.tab//g' )
awk -v name=$name '{OFS="\t"}{print name,$0}' $inTab >> $outDir"all.rRNAcount.tab"

done

#aggregate rRNA counts into one file.

for inTab in $outCountDirAllGenes*.tab ;

do
name=$(basename $inTab | sed 's/.allGenescount.tab//g' )
awk -v name=$name '{OFS="\t"}{print name,$0}' $inTab >> $outDir"all.allGenescount.tab"

done
