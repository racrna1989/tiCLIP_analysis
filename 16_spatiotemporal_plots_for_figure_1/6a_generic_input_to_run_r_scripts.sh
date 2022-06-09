#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 80G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio

#conda install r-base==3.5.1 --force-reinstall --yes

WD="/home/racrna/4_xiCLIP"
cd $WD

OUTDIR="6_mRNACoverage_spatiotemporalMaps/figure_1_heatmaps/"
mkdir -p $OUTDIR

SCRIPT="master_scripts/16_spatiotemporal_plots_for_figure_1/process_counts_over_1kb_binned_RNAs_nolimits_increased_size_range.R"


#filepaths for intronic counts 

INFILE1="6_mRNACoverage_spatiotemporalMaps/all_exon_and_intron.hg38_HeLa_Soren.1kbins.200409.sense.counts" 
OUTFILE1=$OUTDIR"all_exon_and_intron.hg38_HeLa_Soren.1kbins.200409.processed_without_limits_sizeRange_more_sizeRanges.200907.tab"



Rscript $SCRIPT $INFILE1 $OUTFILE1









