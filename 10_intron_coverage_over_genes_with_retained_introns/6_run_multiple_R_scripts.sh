#!/bin/bash

WD="/home/racrna/4_xiCLIP/"
cd $WD
INDIR="master_scripts/10_intron_coverage_over_genes_with_retained_introns/"
INPUTGENERICSCRIPT=$INDIR"5_generic_input_to_run_r_scripts.sh"
OUTDIR=$INDIR"4_induv_Rscripts/"

mkdir -p $OUTDIR

for DESC in ${INDIR}calculate*R ; do 

NAME=$(basename $DESC | sed 's/.sh//g' )

JOBFILE=$OUTDIR"5_"$NAME".sh"

sed "s|CHANGEME|$DESC|g" $INPUTGENERICSCRIPT > $JOBFILE

fBatch $JOBFILE

done
