#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

WD="/home/racrna/4_xiCLIP/"
cd $WD

GENOME="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.chr.sizes"
RPMFACTOR="4_geneCount/allCounts.txt"


for BED_FILEPATH in 5_bedGraphs_and_derivatives/*/*.bed
do

	#This is to construct the output name
	OUTPUT=$(basename $BED_FILEPATH  |sed 's/.bed//g')
	#this is to wrangle sample name to enable selection of scaling factor from RPMFACTOR list
	SAMPLE=$(basename $OUTPUT | sed 's/\./\t/g' | cut -f 1)
	#grep SAMPLE from file, but first turn tab del file into ; del file, remove val. Select Total read number
	TOTALREADS=$(grep $SAMPLE <(sed 's/\t/:/g' $RPMFACTOR | grep -v val) | sed 's/:/\t/g' | cut -f 2) 
	MILL=1000000
	#divide by a million, and save scaling to variable
	SCALE=$(echo "scale=3 ; $MILL / $TOTALREADS" | bc | awk '{printf "%f", $0}' )
	

	#make bedgraph file NOTE grep -v K and G removes random chromosomes 
	bedtools genomecov -bg -split -scale $SCALE -strand + -i $BED_FILEPATH -g $GENOME | grep -v "G" | grep -v "K" > ${BED_FILEPATH/.bed/.RPMScaled.fwd.bedgraph} 
	bedtools genomecov -bg -split -scale $SCALE -strand - -i $BED_FILEPATH -g $GENOME | grep -v "G" | grep -v "K" > ${BED_FILEPATH/.bed/.RPMScaled.rev.bedgraph}
	wait

	bedGraphToBigWig ${BED_FILEPATH/.bed/.RPMScaled.fwd.bedgraph} $GENOME ${BED_FILEPATH/.bed/.RPMScaled.fwd.bw} 
	bedGraphToBigWig ${BED_FILEPATH/.bed/.RPMScaled.rev.bedgraph} $GENOME ${BED_FILEPATH/.bed/.RPMScaled.rev.bw}
	wait

done
