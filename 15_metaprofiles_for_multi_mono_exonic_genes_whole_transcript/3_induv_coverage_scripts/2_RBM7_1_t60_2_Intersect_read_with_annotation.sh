#!/bin/bash
#SBATCH --partition express
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

bedAnnotation="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed"

outDir="8_rna_binding_maps/7_relative_distance_to_TSS_whole_transcript/"
mkdir -p $outDir

inBedDir="5_bedGraphs_and_derivatives/5primepos/RBM7_1_t60"


bedCoverage.relDistTSS.sh $bedAnnotation $inBedDir $outDir

echo "completed bedCoverage"

