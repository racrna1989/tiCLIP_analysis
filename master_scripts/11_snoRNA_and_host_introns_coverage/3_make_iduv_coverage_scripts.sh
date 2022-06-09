#!/bin/bash

WD="/home/racrna/4_xiCLIP/"
cd $WD

OUTDIR="master_scripts/11_snoRNA_and_host_introns_coverage/3_induv_coverage_scripts/"

mkdir -p $OUTDIR

for f in 5_bedGraphs_and_derivatives/read1/*.rRNAScaled.plus.bedgraph ; do 

NAME=$(basename $f | sed 's/\./\t/g' | cut -f 1 )

DESC=$(echo $NAME"*rRNAScaled") 

JOBFILE=$OUTDIR"2_"$NAME"intersect_reads_with_intron_exon_regions.sh"

sed "s/CHANGEME/$DESC/g" master_scripts/11_snoRNA_and_host_introns_coverage/2_intersect_reads_with_intron_exon_regions.sh > $JOBFILE

fBatch $JOBFILE

done