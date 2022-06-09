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
WD="/home/racrna/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/2_Wilson_EJC_iCLIP/"
cd $WD

#annotation files 

ANNO="/home/racrna/4_xiCLIP/annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS.bed"
INBEDDIR="6_bed_files_and_derivatives/5primepos/"
OUTDIR="9_relDistToTSS_mature_210221/1_intersections/"

mkdir -p ${INBEDDIR}

mv 6_bed_files_and_derivatives/*5primepos* ${INBEDDIR}

#functions 

bedCoverage.relDistTSS.sh $ANNO $INBEDDIR $OUTDIR

echo "completed bedCoverage"

outName=$(basename $ANNO | sed 's/bed/count/g')

cat ${OUTDIR}* > ${OUTDIR}xiCLIP_${outName}
