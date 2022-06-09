#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio

WD="/home/racrna/4_xiCLIP/6_mRNACoverage_spatiotemporalMaps/200713-200810-intronic-mRNA-2-1kb/"

cd $WD
pwd

Rscript /home/racrna/4_xiCLIP/6_mRNACoverage_spatiotemporalMaps/200713-200810-intronic-mRNA-2-1kb/200721-intronic-calculate-speed-5.R


