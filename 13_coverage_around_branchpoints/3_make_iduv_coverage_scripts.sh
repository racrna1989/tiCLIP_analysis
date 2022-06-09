#!/bin/bash

WD="/home/racrna/4_xiCLIP/"
cd $WD

OUTDIR="master_scripts/13_coverage_around_branchpoints/3_induv_coverage_scripts/"
INSCRIPT="master_scripts/13_coverage_around_branchpoints/2_intersect_reads_with_branchpoint_regions.sh"
mkdir -p $OUTDIR

for f in 5_bedGraphs_and_derivatives/read1/*RBM7*.rRNAScaled.plus.bedgraph ; do 

NAME=$(basename $f | sed 's/\./\t/g' | cut -f 1 )

DESC=$(echo $NAME"*rRNAScaled") 

OUTSUFFIX=$(basename $INSCRIPT)

JOBFILE=$OUTDIR"2_"$NAME"_"$OUTSUFFIX

sed "s/CHANGEME/$DESC/g" $INSCRIPT > $JOBFILE

fBatch $JOBFILE

done