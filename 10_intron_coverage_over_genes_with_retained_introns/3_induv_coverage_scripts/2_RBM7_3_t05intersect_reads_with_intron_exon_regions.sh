#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP


WD="/home/racrna/4_xiCLIP/"
cd $WD


OUTDIR="8_rna_binding_maps/2_intron_exon_junction_coverage_retained_introns/"

mkdir -p $OUTDIR

INDIR="5_bedGraphs_and_derivatives/"

inAnnotation5end="annotationFiles/hg38_HeLa_transcriptsContainingRI.annotated.v4.200228.5end.100ntupdown.binN201.bed"
inAnnotation3end="annotationFiles/hg38_HeLa_transcriptsContainingRI.annotated.v4.200228.3end.100ntupdown.binN201.bed"

#Read2
outDirRead2=$OUTDIR"read2-bg-rRNAscaled-200228/"
inBedGraphRead2=$INDIR"read2/RBM7_3_t05*rRNAScaled"
mkdir -p $outDirRead2 
bedGraphCoverage.v5.sh $inAnnotation5end $inBedGraphRead2 $outDirRead2 &
bedGraphCoverage.v5.sh $inAnnotation3end $inBedGraphRead2 $outDirRead2 &

#5primepos
outDir5primepos=$OUTDIR"5primepos-bg-rRNAscaled-200228/"
inBedGraph5primepos=$INDIR"5primepos/RBM7_3_t05*rRNAScaled"
mkdir -p $outDir5primepos 
bedGraphCoverage.v5.sh $inAnnotation5end $inBedGraph5primepos $outDir5primepos &
bedGraphCoverage.v5.sh $inAnnotation3end $inBedGraph5primepos $outDir5primepos &

#3endOfRead2
outDir3endOfRead2=$OUTDIR"3endOfRead2-bg-rRNAscaled-200228/"
inBedGraph3endOfRead2=$INDIR"3endOfRead2/RBM7_3_t05*rRNAScaled"
mkdir -p $outDir3endOfRead2
bedGraphCoverage.v5.sh $inAnnotation5end $inBedGraph3endOfRead2 $outDir3endOfRead2 &
bedGraphCoverage.v5.sh $inAnnotation3end $inBedGraph3endOfRead2 $outDir3endOfRead2 &

#Read1
outDirRead1=$OUTDIR"read1-bg-rRNAscaled-200228/"
inBedGraphRead1=$INDIR"read1/RBM7_3_t05*rRNAScaled"
mkdir -p $outDirRead1
bedGraphCoverage.v5.sh $inAnnotation5end $inBedGraphRead1 $outDirRead1 &
bedGraphCoverage.v5.sh $inAnnotation3end $inBedGraphRead1 $outDirRead1 &

wait