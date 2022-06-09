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
rRNAFACTOR="4_geneCount/rRNAFactor.tab"


for BED_FILEPATH in 5_bedGraphs_and_derivatives/*/*.bed

do

	#This is to construct the output name
	OUTPUT=$(basename $INPUT  |sed 's/.bed//g')
	#this is to wrangle sample name to enable selection of scaling factor from rRNAFACTOR list
	SAMPLE=$(basename $INPUT | sed 's/\./\t/g' | cut -f 1)
	SCALE=$(grep $SAMPLE <(cat $rRNAFACTOR | grep -v val) | sed 's/:/\t/g' | cut -f 2) 





	#make bedgraph file NOTE grep -v K and G removes random chromosomes 
	bedtools genomecov -bg -split -scale $SCALE -strand + -i $BED_FILEPATH -g $GENOME | grep -v "G" | grep -v "K" > ${BED_FILEPATH/.bed/.rRNAScaled.fwd.bedgraph} 
	bedtools genomecov -bg -split -scale $SCALE -strand - -i $BED_FILEPATH -g $GENOME | grep -v "G" | grep -v "K" > ${BED_FILEPATH/.bed/.rRNAScaled.rev.bedgraph}
	wait

	bedGraphToBigWig ${BED_FILEPATH/.bed/.rRNAScaled.fwd.bedgraph} $GENOME ${BED_FILEPATH/.bed/.rRNAScaled.fwd.bw} 
	bedGraphToBigWig ${BED_FILEPATH/.bed/.rRNAScaled.rev.bedgraph} $GENOME ${BED_FILEPATH/.bed/.rRNAScaled.rev.bw}
	wait

done
