#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 00:02:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio

WD1="/home/racrna/4_xiCLIP/master_scripts/4_count_over_with_genes_and_total_reads/"
WD2="/home/racrna/4_xiCLIP/"

cd $WD2

Rscript $WD1"process_rRNA_factors.R" 4_geneCount/all.rRNAcount.tab 4_geneCount/rRNAFactor.tab

sed -i 's/ /:/g' 4_geneCount/rRNAFactor.tab
