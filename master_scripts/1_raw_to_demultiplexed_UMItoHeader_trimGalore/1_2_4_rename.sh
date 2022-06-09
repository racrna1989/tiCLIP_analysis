#!/bin/bash 

#rename files to have read number before prefix of .fastq. This is required for trim_galore

WD="/home/racrna/4_xiCLIP/2_demultiplexing/3_RT_L3_demultiplexed/"
cd $WD
for f in *R1.processed*; do rename .fastq _1.fastq $f; done
for f in *R2.processed*; do rename .fastq _2.fastq $f; done
