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

SCRIPT="CHANGEME"


ANNOBED="annotationFiles/hg38_HeLa_transcriptsContainingRI.annotated.v4.200228.bed"
EXPRVEC="annotationFiles/log2_mean_cov_RNAseq_TTseq.RData"

for inCount in 8_rna_binding_maps/2_intron_exon_junction_coverage_retained_introns/*200228/*3end.100nt*counts; do 

	#wrangle input filepaths for 3 and 5 end
	THREE_END=$inCount
	FIVE_END=$(echo $inCount | sed 's/\.3end\.100nt/\.5end\.100nt/g')
	
	#get file descriptor from input annotation
	FILEDESC=$(basename $SCRIPT | sed 's/.R//g ')
	
	#edit input file path to make output file path for output file
	OUT_NAME_FILE=$(basename $inCount | sed 's/\./\t/g' | awk -v FILEDESC=$FILEDESC '{OFS="\t"}{print $1"_"$2"_"$3"."FILEDESC".tab"}')
	OUT_NAME_DIR=$(dirname $inCount)
	OUT_NAME="$OUT_NAME_DIR/$OUT_NAME_FILE"
	
	echo $THREE_END 
	echo $FIVE_END
	
	Rscript $SCRIPT $ANNOBED $EXPRVEC $THREE_END $FIVE_END $OUT_NAME

done
