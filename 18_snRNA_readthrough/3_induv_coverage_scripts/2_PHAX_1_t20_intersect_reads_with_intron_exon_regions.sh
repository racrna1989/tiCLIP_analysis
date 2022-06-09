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

OUTDIR="9_snRNA_readthrough/2_snRNA_coverage_rel_dist_3end/"

mkdir -p $OUTDIR

INDIR="5_bedGraphs_and_derivatives/"

inAnnotation="annotationFiles/snRNA.GRCh38101_and_Soren_anno.interval100ntup2500down.bed"

#Read2
outDirRead2=${OUTDIR}"read2-bg-rRNAscaled-200228/"
inBedGraphRead2=$INDIR"read2/PHAX_1_t20"
mkdir -p $outDirRead2 
bedCoverage.relDistTSS.sh $inAnnotation $inBedGraphRead2 $outDirRead2 &

#5primepos
outDir5primepos=$OUTDIR"5primepos-bg-rRNAscaled-200228/"
inBedGraph5primepos=$INDIR"5primepos/PHAX_1_t20"
mkdir -p $outDir5primepos 
bedCoverage.relDistTSS.sh $inAnnotation $inBedGraph5primepos $outDir5primepos &

#3endOfRead2
outDir3endOfRead2=$OUTDIR"3endOfRead2-bg-rRNAscaled-200228/"
inBedGraph3endOfRead2=$INDIR"3endOfRead2/PHAX_1_t20"
mkdir -p $outDir3endOfRead2
bedCoverage.relDistTSS.sh $inAnnotation $inBedGraph3endOfRead2 $outDir3endOfRead2 &

#Read1
outDirRead1=$OUTDIR"read1-bg-rRNAscaled-200228/"
inBedGraphRead1=$INDIR"read1/PHAX_1_t20"
mkdir -p $outDirRead1
bedCoverage.relDistTSS.sh $inAnnotation $inBedGraphRead1 $outDirRead1 &

wait