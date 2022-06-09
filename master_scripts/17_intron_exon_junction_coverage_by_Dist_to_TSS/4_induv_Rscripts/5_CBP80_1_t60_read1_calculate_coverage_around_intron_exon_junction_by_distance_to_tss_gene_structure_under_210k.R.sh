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

#conda install r-base==3.5.1 --force-reinstall --yes

WD="/home/racrna/4_xiCLIP/"
cd $WD

SCRIPT="calculate_coverage_around_intron_exon_junction_by_distance_to_tss_gene_structure_under_210k.R"


ANNOBED="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS.bed"

EXPRVEC="annotationFiles/log2_mean_cov_RNAseq_TTseq.RData"

for inCount in 8_rna_binding_maps/4_intron_exon_junction_coverage_by_exon_context/*200228/CBP80_1_t60_read1*3end.100nt*counts; do 

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

