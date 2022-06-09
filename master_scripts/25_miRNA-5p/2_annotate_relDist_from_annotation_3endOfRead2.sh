#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/"

cd ${WD}



ANNO="annotationFiles/miRNA-5p_annotated.5end.100ntupdown.binN201.bed"

#3 calculate rel distance

INBEDGRAPH="5_bedGraphs_and_derivatives/3endOfRead2/"
OUTDIR="8_rna_binding_maps/7_miRNA-5p_5end/2_coverage_over_miRNA-5p/3endOfRead2/"
mkdir -p ${OUTDIR}
	
bedGraphCoverage.v1.rRNAScaled.sh ${ANNO} ${INBEDGRAPH} ${OUTDIR}