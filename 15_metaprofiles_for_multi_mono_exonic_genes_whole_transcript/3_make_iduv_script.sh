#!/bin/bash

WD="/home/racrna/4_xiCLIP/"
cd $WD

OUTDIR="master_scripts/15_metaprofiles_for_multi_mono_exonic_genes_whole_transcript/3_induv_coverage_scripts/"
INSCRIPT="master_scripts/15_metaprofiles_for_multi_mono_exonic_genes_whole_transcript/2_Intersect_read_with_annotation.sh"
mkdir -p $OUTDIR

for f in 5_bedGraphs_and_derivatives/5primepos/*.bed ; do 

NAME=$(basename $f | sed 's/\./\t/g' | cut -f 1 )


OUTSUFFIX=$(basename $INSCRIPT)

JOBFILE=$OUTDIR"2_"$NAME"_"$OUTSUFFIX

sed "s/CHANGEME/$NAME/g" $INSCRIPT > $JOBFILE

fBatch $JOBFILE

done