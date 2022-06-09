#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/"

cd ${WD}


OUTDIR="8_rna_binding_maps/7_miRNA-5p_5end/2_coverage_over_miRNA-5p/"
cd ${OUTDIR}

for f in */*.counts
do 
READTYPE=$(basename $f | cut -d "." -f 2)
OUTFILE=$(basename $f | cut -d "." -f 4)

awk -v READTYPE=$READTYPE -v OUTFILENAME=$OUTFILENAME '{OFS="\t"}{
$1=$1"_"READTYPE; print $0
}' $f >> "xiCLIP_"${READTYPE}"_"${OUTFILE}".tab"
done