#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

# This file will preform counts over all the datasets. Using specifically read2 with the strand stop, 
#taken from /home/racrna/3_xiCLIP/16_geneCounts_intronExonRatio/

WD="/home/racrna/4_xiCLIP/"

cd $WD

ANNOBED="/home/racrna/faststorage/GRCh38/HeLa/custom/hg38_HeLa_Soren.1kbins.200409.bed"
RRNAFACTOR="4_geneCount/rRNAFactor.tab"

OUTDIR="4_geneCount/binnedGeneCounts-1kb/"
INDIR="5_bedGraphs_and_derivatives/read2/"

mkdir -p $OUTDIR

#count overlap
for INBED in ${INDIR}*read2.bed 
do
	OUTNAME=$(basename ${INBED} | sed 's/.bed//g')
	bedtools intersect -c -s -a ${ANNOBED} -b ${INBED} > ${OUTDIR}${OUTNAME}".hg38HeLaSoren1kbins.200409.counts"
	
done

#normalise to rRNA scaling factor

for f in ${OUTDIR}*.hg38HeLaSoren1kbins.200409.counts; 
do 
	OUTNAME=$(basename $f | sed 's/.read2.hg38HeLaSoren1kbins.200409.counts//g'); 
	echo $OUTNAME ; 
	SCALE=$(grep $OUTNAME $RRNAFACTOR | awk '{print $2}') ; 
	awk -v sF=$SCALE -v name=$OUTNAME '{OFS="\t"}{$NF=($NF*sF);print name,$0}' $f > ${OUTDIR}${OUTNAME}"_hg38HeLaSoren1kbins.NormTorRNAFactor.counts"
done
