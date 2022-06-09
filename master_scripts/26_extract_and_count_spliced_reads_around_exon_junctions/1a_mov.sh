#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

#set wd
WD="/home/racrna/4_xiCLIP/"
cd ${WD}

MASTEROUTDIR="12_extract_splice_and_non_spliced_reads_close_to_exon_junctions/"
OUTDIR=${MASTEROUTDIR}"1_spliced_and_not_spliced_bams/"

mv *.bam ${OUTDIR}