#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/annotationFiles"
cd $WD

cat ~/3_xiCLIP/13_heatmaps_MetaGeneAnalysis/annotationFiles/tmp/BP-annotated-metadata.hg38.100ntupdown.binN201_* | sort -k1,1 -k2,2n > BP-annotated-metadata.hg38.100ntupdown.binN201.bed
