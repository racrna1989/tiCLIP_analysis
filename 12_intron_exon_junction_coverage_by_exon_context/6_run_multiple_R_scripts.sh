#!/bin/bash

WD="/home/racrna/4_xiCLIP/"
cd $WD

INDIR="master_scripts/12_intron_exon_junction_coverage_by_exon_context/"
INPUTGENERICSCRIPT=$INDIR"5_generic_input_to_run_r_scripts.sh"
OUTDIR=$INDIR"4_induv_Rscripts/"
ANNOTATION_FILE_PATH="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS.bed"


mkdir -p $OUTDIR

for DESC in ${INDIR}calculate*R ; do 

NAME=$(basename $DESC | sed 's/.sh//g' )

JOBFILE=$OUTDIR"5_"$NAME".sh"

sed "s|CHANGEME|$DESC|g;s|ANNOTATIONFILEPATH|$ANNOTATION_FILE_PATH|g" $INPUTGENERICSCRIPT > $JOBFILE

fBatch $JOBFILE

done



