#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio


WD="/home/racrna/4_xiCLIP/"
cd $WD
ANNOBED="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.bed"
SCRIPT="master_scripts/9_all_introns_coverage_plots/4_wrangle_intron_counts_and_norm_to_gene_expression_and_number_of_introns.R"
EXPRVEC="annotationFiles/log2_mean_cov_RNAseq_TTseq.RData"
for inCount in 8_rna_binding_maps/1_intron_exon_junction_coverage/*/*3end.100nt*; do 

THREE_END=$inCount

FIVE_END=$(echo $inCount | sed 's/\.3end\.100nt/\.5end\.100nt/g')

OUT_NAME=$(echo $inCount | sed 's/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.3end.100ntupdown.binN201.sense.counts/intron_coverage_100ntupdown_sense.tab/g;s/.rRNAScaled/_rRNAScaled/g')

echo $THREE_END 

echo $FIVE_END

Rscript $SCRIPT $ANNOBED $EXPRVEC $THREE_END $FIVE_END $OUT_NAME

done