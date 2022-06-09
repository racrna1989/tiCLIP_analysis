#!/bin/bash
#SBATCH --partition express
#SBATCH --mem-per-cpu 40G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#set WD
WD="/home/racrna/4_xiCLIP/"
cd ${WD}

#annotation files 

ANNO="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS.bed"
INBEDDIR="5_bedGraphs_and_derivatives/5primepos/"
OUTDIR="13_relDistToTSS_mature_200416/1_intersections/"


#functions 

bedCoverage.relDistTSS.sh $ANNO $INBEDDIR $OUTDIR

echo "completed bedCoverage"

outName=$(basename $ANNO | sed 's/bed/count/g')

cat ${OUTDIR}* > ${OUTDIR}xiCLIP_${outName}



