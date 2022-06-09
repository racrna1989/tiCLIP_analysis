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

ANNOBED="annotationFiles/snoRNAs.GRCh38andrefGene.mature.bed"
SCRIPT="master_scripts/11_snoRNA_and_host_introns_coverage/4_wrangle_retained_intron_coverage_files.R"
EXPRVEC="annotationFiles/log2_mean_cov_RNAseq_TTseq.RData"
HOST_GENE_INDEX_FILEPATH="annotationFiles/snoRNA_containing_transcripts_index.tab"

for inCount in 8_rna_binding_maps/3_snoRNA_and_host_intron_coverage_at_ends_of_annotated_genes/*/*counts; do 

IN_FILEPATH=$inCount

OUT_NAME=$(echo $inCount | sed 's/\./\t/g' | awk '{OFS="\t"}{print $1"_"$2"_"$3"."$4".tab"}' )

echo $IN_FILEPATH 

Rscript $SCRIPT $ANNOBED $EXPRVEC $IN_FILEPATH $HOST_GENE_INDEX_FILEPATH $OUT_NAME

done



