#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk


source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

WD="/home/racrna/4_xiCLIP/annotationFiles"
cd $WD

ANNO="hg38_HeLa_transcriptsContainingRI.annotated.v4.200228.bed"

for FILE in $ANNO; do

OUTFILE=$(echo $FILE | sed 's/.bed//g;s/\./_/g')
get5primeFlank.bash $ANNO 100 100 > $OUTFILE".100nt_5prime.bed"
get3primeFlank.bash $ANNO 100 100 > $OUTFILE".100nt_3prime.bed"
getWindows_3end.v3.sh $ANNO 100 100
getWindows_5end.v3.sh $ANNO 100 100
done

