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


OUTDIR="8_rna_binding_maps/3_snoRNA_and_host_intron_coverage_at_ends_of_annotated_genes/"

mkdir -p $OUTDIR

INDIR="5_bedGraphs_and_derivatives/"

inAnnotation="annotationFiles/hg38_snoRNAs_and_host_introns.100ntupdown.binN201.bed"

#Read2
outDirRead2=$OUTDIR"read2-bg-rRNAscaled-200228/"
inBedGraphRead2=$INDIR"read2/PHAX_1_t20*rRNAScaled"
mkdir -p $outDirRead2 
bedGraphCoverage.v5.sh $inAnnotation $inBedGraphRead2 $outDirRead2 &

#5primepos
outDir5primepos=$OUTDIR"5primepos-bg-rRNAscaled-200228/"
inBedGraph5primepos=$INDIR"5primepos/PHAX_1_t20*rRNAScaled"
mkdir -p $outDir5primepos 
bedGraphCoverage.v5.sh $inAnnotation $inBedGraph5primepos $outDir5primepos &

#3endOfRead2
outDir3endOfRead2=$OUTDIR"3endOfRead2-bg-rRNAscaled-200228/"
inBedGraph3endOfRead2=$INDIR"3endOfRead2/PHAX_1_t20*rRNAScaled"
mkdir -p $outDir3endOfRead2
bedGraphCoverage.v5.sh $inAnnotation $inBedGraph3endOfRead2 $outDir3endOfRead2 &

#Read1
outDirRead1=$OUTDIR"read1-bg-rRNAscaled-200228/"
inBedGraphRead1=$INDIR"read1/PHAX_1_t20*rRNAScaled"
mkdir -p $outDirRead1
bedGraphCoverage.v5.sh $inAnnotation $inBedGraphRead1 $outDirRead1 &

wait