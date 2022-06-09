#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/"

cd $WD

#bedGraphCoverage.specialisedForSpatioTemporal_intron_counts.v3.sh annotationFiles/hg38_HeLa_Soren.1kbins.200409.bed 5_bedGraphs_and_derivatives/5primepos/ 6_mRNACoverage_spatiotemporalMaps/200713-200810-intronic-mRNA-2-1kb/ &
#bedGraphCoverage.specialisedForSpatioTemporal_intron_counts.v3.sh annotationFiles/hg38_HeLa_Soren.10kbins.200409.bed 5_bedGraphs_and_derivatives/5primepos/ 6_mRNACoverage_spatiotemporalMaps/200713-200810-intronic-mRNA-2-10kb/ &

bedGraphCoverage.specialisedForSpatioTemporal_exon_counts.v3.sh annotationFiles/hg38_HeLa_Soren.1kbins.200409.bed 5_bedGraphs_and_derivatives/5primepos/ 6_mRNACoverage_spatiotemporalMaps/200713-200810-exonic-mRNA-2-1kb/ 
