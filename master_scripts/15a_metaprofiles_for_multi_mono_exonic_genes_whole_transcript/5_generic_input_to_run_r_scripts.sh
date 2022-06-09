#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 40G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio


WD="/home/racrna/4_xiCLIP/"
cd $WD

INDIR="8_rna_binding_maps/8_relative_distance_to_TSS_whole_transcript_increased_range"
ANNOTATION_BED_FILEPATH="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed"
SCRIPT="master_scripts/15a_metaprofiles_for_multi_mono_exonic_genes_whole_transcript/4_normalise_and_scale_relative_dist_from_TSS_whole_transcript.R"


rRNA_FACTOR_FILEPATH="4_geneCount/allCounts.txt"
RPM_FACTOR_FILEPATH="4_geneCount/rRNAFactor.tab"

EXPRVEC="annotationFiles/log2_mean_cov_RNAseq_TTseq.RData"

for DATA_FRAME_FILEPATH in $INDIR/*counts; do 
	
	#get file descriptor from input annotation
	FILEDESC=$(basename $SCRIPT | sed 's/.R//g ')
	
	#edit input file path to make output file path for output file
	OUT_NAME_FILE=$(basename $DATA_FRAME_FILEPATH | sed 's/_/\t/g' | awk -v FILEDESC=$FILEDESC '{OFS="\t"}{print $1"_"$2"_"$3"."FILEDESC".tab"}')
	OUT_NAME_DIR=$(dirname $DATA_FRAME_FILEPATH)
	OUTFILE_PATH="$OUT_NAME_DIR/$OUT_NAME_FILE"
	
	echo "processing $OUT_NAME_FILE"
	
	if [ -s "$OUTFILE_PATH" ] ; 
	then
		echo "file exists and is not empty, skipping"
	else
		Rscript $SCRIPT $DATA_FRAME_FILEPATH $ANNOTATION_BED_FILEPATH $rRNA_FACTOR_FILEPATH $RPM_FACTOR_FILEPATH $OUTFILE_PATH
	fi
	
done

#aggregate files into one long table


for DATA_FRAME_FILEPATH in $INDIR/*tab; do 

#get file descriptor from input annotation
	FILEDESC=$(basename $SCRIPT | sed 's/.R//g ')
	
	#edit input file path to make output file path for output file
	OUT_NAME_FILE="xiCLIP_all."$FILEDESC".count"
	OUT_NAME_DIR=$(dirname $DATA_FRAME_FILEPATH)
	OUTFILE_PATH="$OUT_NAME_DIR/$OUT_NAME_FILE"
	
	tail -n +2 $DATA_FRAME_FILEPATH >> $OUTFILE_PATH
	
done





