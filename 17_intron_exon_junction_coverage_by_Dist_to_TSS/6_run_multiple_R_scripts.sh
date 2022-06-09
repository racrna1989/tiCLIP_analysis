#!/bin/bash

WD="/home/racrna/4_xiCLIP/"
cd $WD

#input the filepath from 4_xiCLIP to the dir containing scripts for this specific operation

INDIR="master_scripts/17_intron_exon_junction_coverage_by_Dist_to_TSS/"
INPUTGENERICSCRIPT=$INDIR"5_generic_input_to_run_r_scripts.sh"
OUTDIR=$INDIR"4_induv_Rscripts/"
ANNOTATION_FILE_PATH="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS.bed"


mkdir -p $OUTDIR

for DESC in ${INDIR}calculate*R ; do 
for SAMPLE_NAME_FILEPATH in 8_rna_binding_maps/4_intron_exon_junction_coverage_by_exon_context/*200228/*3end.100nt*counts; do

SCRIPT=$(basename $DESC | sed 's/.sh//g' )

SAMPLE_NAME=$(basename $SAMPLE_NAME_FILEPATH | sed 's/\./\t/g' | awk '{OFS="\t"}{print $1"_"$2}')

JOBFILE=$OUTDIR"5_"$SAMPLE_NAME"_"$SCRIPT".sh"

sed "s|CHANGEME|$SCRIPT|g;s|ANNOTATIONFILEPATH|$ANNOTATION_FILE_PATH|g;s|SAMPLE|$SAMPLE_NAME|g" $INPUTGENERICSCRIPT > $JOBFILE
echo $JOBFILE
#fBatch $JOBFILE

done
done



