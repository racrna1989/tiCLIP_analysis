#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/"

cd ${WD}

OUTDIR="8_rna_binding_maps/7_miRNA-5p_5end/"

mkdir -p ${OUTDIR}

ANNO1="/home/racrna/faststorage/GRCh38/miRNA/miRNA-5p_annotated.5end.100ntupdown.binN201.bed"
ANNO2="/home/racrna/faststorage/GRCh38/miRNA/miRNA_5p_5end.bed"

cd annotationFiles/

#ln -s ${ANNO1} . 
ln -s ${ANNO2} . 
