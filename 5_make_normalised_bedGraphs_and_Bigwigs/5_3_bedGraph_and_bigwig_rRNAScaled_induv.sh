#!/bin/bash

#bash file to echo a script for each induvidual bed file to make bedgraphs and bigwigs that are scaled to rRNA counts.

WD="/home/racrna/4_xiCLIP/"
cd $WD

jobDir="jobs/5_make_normalised_bedGraphs_and_Bigwigs/"
rRNAFACTOR="4_geneCount/rRNAFactor.tab"
mkdir -p $jobDir

for INPUT in 5_bedGraphs_and_derivatives/*/*.bed; do 
	#This is to construct the output name
	OUTPUT=$(basename $INPUT  |sed 's/.bed//g')
	#this is to wrangle sample name to enable selection of scaling factor from rRNAFACTOR list
	SAMPLE=$(basename $INPUT | sed 's/\./\t/g' | cut -f 1)
	SCALE=$(grep $SAMPLE <(cat $rRNAFACTOR | grep -v val) | sed 's/:/\t/g' | cut -f 2) 

	#create an empty file in the jobs directory to save the mapping script for later reference.
	jobFile="${jobDir}5_3_bedGraph_and_bigwig_rRNAScaled__dz_option${OUTPUT}.job"

## here we do a four part echo to echo the script and the variable into the jobFile created above. 
#I used three seperate echos to allow the echoing of the " in each segment 
#FIRST


echo '#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

WD="/home/racrna/4_xiCLIP/"
cd $WD

GENOME="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.chr.sizes"
' > $jobFile

#SECOND
echo 'for BED_FILEPATH in '$INPUT >> $jobFile

#THIRD
echo '
do

SCALE='${SCALE} >> $jobFile

#FORTH - new piece is the dz option on the genomecove and the additional awk script at the end of the pipe. 
echo '
#make bedgraph file NOTE grep -v K and G removes random chromosomes 
bedtools genomecov -dz -split -scale $SCALE -strand + -i $BED_FILEPATH -g $GENOME | grep -v "G" | grep -v "K" | awk \'{OFS="\t"}{print $1,$2,$2+1,$3}\' > ${BED_FILEPATH/.bed/.rRNAScaled.dzOption.fwd.bedgraph} & 
bedtools genomecov -dz -split -scale $SCALE -strand - -i $BED_FILEPATH -g $GENOME | grep -v "G" | grep -v "K" | awk \'{OFS="\t"}{print $1,$2,$2+1,$3}\' > ${BED_FILEPATH/.bed/.rRNAScaled.dzOption.rev.bedgraph} &

wait

bedGraphToBigWig ${BED_FILEPATH/.bed/.rRNAScaled.dzOption.fwd.bedgraph} $GENOME ${BED_FILEPATH/.bed/.rRNAScaled.dzOption.fwd.bw} &
bedGraphToBigWig ${BED_FILEPATH/.bed/.rRNAScaled.dzOption.rev.bedgraph} $GENOME ${BED_FILEPATH/.bed/.rRNAScaled.dzOption.rev.bw} &

wait

done' >> $jobFile

##the job file is then submitted using the fBatch wrapper. The output file for this should contain the jobis in a list
fBatch $jobFile

sleep 2

done

