#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/"

cd $WD

ANNOBED="annotationFiles/snRNA_GRCh38101_and_Soren_anno_2500ntextention.bed"
INDIR="5_bedGraphs_and_derivatives/"
OUTDIR="9_snRNA_readthrough/1_binding_densities/"

mkdir -p $OUTDIR

for INBED in ${INDIR}*/*.bed 
do
	
	TMPFILE=${OUTDIR}".${RANDOM}.bed"
	
	touch ${TMPFILE} 
	
	OUTNAME=$(basename ${INBED} | sed 's/.bed//g' | sed 's/\./_/g')
	
	DESC=$(basename ${ANNOBED} | sed 's/.bed//g')
	
	echo "processing :" ${OUTNAME}
		
	#make tmp file of snRNA overlap
	bedtools intersect -s -a ${ANNOBED} -b ${INBED} > ${TMPFILE}
	
	#count overlap between 2500nt extentions of snRNAs
	
	bedtools intersect -c -s -a <(grep "2500nt" ${ANNOBED}) -b $TMPFILE > ${OUTDIR}${OUTNAME}"."${DESC}".counts"
	
	bedtools intersect -c -s -a <(grep "genebody" ${ANNOBED}) -b <(bedtools intersect -s -v -a $TMPFILE -b <(grep "2500nt" ${ANNOBED})) >> ${OUTDIR}${OUTNAME}"."${DESC}".counts"
	
	rm ${TMPFILE}

	
done