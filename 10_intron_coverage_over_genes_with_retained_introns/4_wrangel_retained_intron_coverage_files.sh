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

ANNOBED="annotationFiles/hg38_HeLa_transcriptsContainingRI.annotated.v4.200228.bed"
SCRIPT="master_scripts/10_intron_coverage_over_genes_with_retained_introns/4_wrangel_retained_intron_coverage_files.R"
EXPRVEC="annotationFiles/log2_mean_cov_RNAseq_TTseq.RData"
for inCount in 8_rna_binding_maps/2_intron_exon_junction_coverage_retained_introns/*/*3end.100nt*; do 

THREE_END=$inCount

FIVE_END=$(echo $inCount | sed 's/\.3end\.100nt/\.5end\.100nt/g')

OUT_NAME=$(echo $inCount | sed 's/hg38_HeLa_transcriptsContainingRI.annotated.v4.200228.3end.100ntupdown.binN201.sense.counts/retained_intron_coverage_100ntupdown_sense.tab/g;s/.rRNAScaled/_rRNAScaled/g')

echo $THREE_END 

echo $FIVE_END

Rscript $SCRIPT $ANNOBED $EXPRVEC $THREE_END $FIVE_END $OUT_NAME

done


