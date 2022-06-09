#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 60G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio


WD="/home/racrna/4_xiCLIP/"
cd $WD

SCRIPT="master_scripts/16_spatiotemporal_plots_for_figure_1/process_counts_over_1kb_binned_RNAs_nolimits.R"



#filepaths for exonic 

INFILE2="6_mRNACoverage_spatiotemporalMaps/200713-200810-exonic-mRNA-2-1kb/all_Exons.hg38_HeLa_Soren.1kbins.200409.sense.counts"
OUTFILE2="6_mRNACoverage_spatiotemporalMaps/200713-200810-exonic-mRNA-2-1kb/all_exons.hg38_HeLa_Soren.1kbins.200409.processed_without_limits_sizeRange.200904.tab"




Rscript $SCRIPT $INFILE2 $OUTFILE2





