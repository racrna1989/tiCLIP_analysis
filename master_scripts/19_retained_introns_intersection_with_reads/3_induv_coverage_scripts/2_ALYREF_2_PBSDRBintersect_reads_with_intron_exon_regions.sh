#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 40G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP


WD="/home/racrna/4_xiCLIP/"
cd $WD

OUTDIR="8_rna_binding_maps/2_intron_exon_junction_coverage_retained_introns/2_read_intersections/"

mkdir -p $OUTDIR

INDIR="5_bedGraphs_and_derivatives/"

#annotation files 
five_annotation="annotationFiles/hg38_HeLa_transcriptsContainingRI_annotated_v4_200228.100nt_5prime.bed"
three_annotation="annotationFiles/hg38_HeLa_transcriptsContainingRI_annotated_v4_200228.100nt_3prime.bed"

#5primepos
outDir5primepos=$OUTDIR"5primepos-bg-rRNAscaled-200228/"
inBedGraph5primepos=$INDIR"5primepos/ALYREF_2_PBSDRB"
mkdir -p $outDir5primepos 
bedCoverage.relDistTSS.sh $five_annotation $inBedGraph5primepos $outDir5primepos &
bedCoverage.relDistTSS.sh $three_annotation $inBedGraph5primepos $outDir5primepos &

#3endOfRead2
outDir3endOfRead2=$OUTDIR"3endOfRead2-bg-rRNAscaled-200228/"
inBedGraph3endOfRead2=$INDIR"3endOfRead2/ALYREF_2_PBSDRB"
mkdir -p $outDir3endOfRead2
bedCoverage.relDistTSS.sh $five_annotation $inBedGraph3endOfRead2 $outDir3endOfRead2 &
bedCoverage.relDistTSS.sh $three_annotation $inBedGraph3endOfRead2 $outDir3endOfRead2 &

wait