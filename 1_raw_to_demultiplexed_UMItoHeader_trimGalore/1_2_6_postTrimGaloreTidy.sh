#!/bin/bash

move reports to their own folders 

WD="/home/racrna/4_xiCLIP/2_demultiplexing/4_RT_L3_demultiplexed_trimGalore/"
cd $WD

fqcDir="fastqc/"
trimRepDir="trimmingReport/"

mkdir -p $fqcDir $trimRepDir

mv *fastqc* $fqcDir
mv *trimming_rep* $trimRepDir
