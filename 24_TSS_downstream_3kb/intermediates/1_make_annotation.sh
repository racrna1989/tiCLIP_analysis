#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

wd="/home/racrna/4_xiCLIP/annotationFiles/"
cd $wd


#making updated snRNA annotation file. 

#step make annotation file for 10 nt up and down from 5SS and BP

#PROMPTSPAXTNEXT contains NA at col 3 for some + annotations. Doesnt matter for this as we will take a range from TSS.


grep :::TtT3::: hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed >> tmp.bed
grep :::HtH5::: hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed >> tmp.bed
grep :::nTtT::: hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed >> tmp.bed
grep :::nHtH::: hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed >> tmp.bed

cat PROMPTPAXTNEXT.bed tmp.bed | sort -k1,1 -k2,2n - > PROMPTs_and_HeLa_PROMPTs.bed

rm tmp.bed

#make TSS file for annotation

awk '{OFS="\t"}{
if ($6 == "+"){
$3=$2+1; print $0}
else if ($6 == "-"){
$2=$3-1; print $0}}' PROMPTs_and_HeLa_PROMPTs.bed > PROMPTs_and_HeLa_PROMPTs.TSS.bed


