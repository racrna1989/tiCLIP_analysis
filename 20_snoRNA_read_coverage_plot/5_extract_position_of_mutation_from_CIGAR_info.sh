#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

##################
#
#Purpose - To extract read1 and ultimately extract the mutational profile of the reads 
#290620
#
##################


WD="/home/racrna/4_xiCLIP/"

cd $WD

OUTDIR="10_snoRNA_read_coverage/4_getMutProfile/"

#########
#PART 2##
#########

#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

INBAMDIR=${OUTDIR}"3_read1-snoRNA-bams/"
MUTATIONDIR=${OUTDIR}"4_read_CIGARinfo/"

#mkdir
mkdir -p ${MUTATIONDIR}

#loop to extract mutational information

for INBAM in ${INBAMDIR}*.read1.bam; 
do 
	#edit filename to save ID
	ID=$(basename ${INBAM} | cut -d "." -f 1)
	echo "analysing ${ID}"
	sam2tsv -R ~/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa ${INBAM} | \
	awk '{OFS="\t"}{
	if ($NF !="M"){
	print $3,$7,$1,$4,$NF"|"$8"|"$5}}' > ${MUTATIONDIR}${ID}".CIGARinfo.tab"
done

#outfile
#$chr $refpos $readID $readPOS $mutCode $readB$refbase 

