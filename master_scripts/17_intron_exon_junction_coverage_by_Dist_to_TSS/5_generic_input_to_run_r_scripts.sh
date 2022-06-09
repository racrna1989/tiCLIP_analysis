#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 100G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio

#conda install r-base==3.5.1 --force-reinstall --yes

WD="/home/racrna/4_xiCLIP/"
cd $WD

SCRIPT="process_exon_intron_coverage_file_2.R"


ANNOBED="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS.bed"

EXPRVEC="annotationFiles/log2_mean_cov_RNAseq_TTseq.RData"

8_rna_binding_maps/4_intron_exon_junction_coverage_by_exon_context/xiCLIP.read*.counts; 
do  
READTYPE=$(basename $COUNTS | cut -d "." -f 2 )

Rscript $SCRIPT $ANNOBED $EXPRVEC $COUNTS $READTYPE

done


