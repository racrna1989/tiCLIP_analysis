#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

wd="/home/racrna/faststorage/GRCh38"
OLDANNO="HeLa/custom/snRNA.2500ntext.bed"
NEWANNO="genome/Homo_sapiens.GRCh38.101.gtf"

#making updated snRNA annotation file. 

#step 1 remove snRNAs from new annotation that are already present in old

grep geneBody $OLDANNO

grep snRNA $NEWANNO | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$18,$7}}' | tr -d '";'


#returns nonoverlapping snRNAs from new annotation
bedtools intersect -s -v -a <(grep snRNA $NEWANNO | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$18,$7}}' | tr -d '";') -b <(grep geneBody $OLDANNO) | awk '{OFS="\t"}{print $1,$2,$3,$2":::"$3":::"$4":::"$5,".",$6}' > tmp1.bed

#makes old file the same format as the new one
grep geneBody $OLDANNO | sed 's/:::/\t/g' | awk '{OFS="\t"}{print $1,$2,$3,$4":::"$5":::"$6":::"$7,$(NF-1),$NF}' > tmp2.bed

#cat and sort files 

cat tmp1.bed tmp2.bed > snRNA.GRCh38101_and_Soren_anno.bed
rm tmp1.bed tmp2.bed
cd snRNA/

awk '{OFS="\t"}{if ($6 =="+"){print $1,$3,$3+2500,$4":::2500nt_downstream",$5,$6} else if ($6 =="-"){print $1, $2-2500,$2,$4":::2500nt_downstream",$5,$6}}' snRNA.GRCh38101_and_Soren_anno.bed > extentions.tmp
awk '{OFS="\t"}{$4=$4":::genebody"; print $0}' snRNA.GRCh38101_and_Soren_anno.bed | cat - extentions.tmp | sort -k1,1 -k2,2n | grep -v "KI27" | grep -v "GL00" > snRNA.GRCh38101_and_Soren_anno.2500ntextention.bed

rm extentions.tmp



#ln file to annotation file in 4_xiCLIP

cd /home/racrna/4_xiCLIP/annotationFiles/

ln -s ~/faststorage/GRCh38/snRNA/snRNA.GRCh38101_and_Soren_anno.2500ntextention.bed .


#make annotation file that contains snRNA exteneded by 2500nt 

cd $wd
cd snRNA 

awk '{OFS="\t"}{if ($6 =="+"){print $1,$2-100,$3+2500,$4":::100ntup2500ntdown",$5,$6} else if ($6 =="-"){print $1,$2-2500,$3+100,$4":::100ntup2500ntdown",$5,$6}}' snRNA.GRCh38101_and_Soren_anno.bed > snRNA.GRCh38101_and_Soren_anno.interval100ntup2500down.bed

cd /home/racrna/4_xiCLIP/annotationFiles/

ln -s ~/faststorage/GRCh38/snRNA/snRNA.GRCh38101_and_Soren_anno.interval100ntup2500down.bed . 

