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
mkdir -p $OUTDIR

#filenames and dirs req


ANNODOWNSTREAMSNO="annotationFiles/snoRNAs.GRCh38andrefGene.maturetodowntreamexon.bed"
ANNOMATURESNO="annotationFiles/snoRNAs.GRCh38andrefGene.mature.bed"


regionBed="annotationFiles/snoRNAs.GRCh38andrefGene.maturetodowntreamexon.bed"
genomeFasta="~/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa"

BAMDIR="3_mapping/2_QC_mapped/"

headerDir=${OUTDIR}"2_samHeaders/"

OUTBAMDIR=${OUTDIR}"3_read1-snoRNA-bams/"

readIDDir=${OUTDIR}"1_readIDList/"

INDIR="10_snoRNA_read_coverage/3_map_snoRNA_to_read/"

#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#make dirs if not already present
mkdir -p $BAMDIR $headerDir $OUTBAMDIR $readIDDir

for INBED in ${INDIR}*.bed; 
do 
#edit filename to save ID
ID=$(basename ${INBED} | cut -d "." -f 1)
echo "analysing ${ID}"
#print out readIDs from bed file, sort and save in readIDlist dir
awk '{print $4}' ${INBED} | sort -k1,1 | uniq > ${readIDDir}${ID}".readIDs.txt"
echo "printing header for read 1 $ID"	
#save header from original read1 bamfile.
samtools view -f 64 -H ${BAMDIR}${ID}".SoUmiDedupRemSec.bam" > $headerDir${ID}".SoUmiDedupRemSec.header" 

#mktmp file
tmp1="$(mktemp racrna.XXXXXX.tmp)"
tmp2="$(mktemp racrna.XXXXXX.tmp)"
echo `ls $tmp1`
echo `ls $tmp2`

echo "grepping IDs from sam"
#extract read1 of reads that overlap bam region file. Grep reads, cat header and extracted reads
samtools view -f 64 -L ${regionBed} ${BAMDIR}${ID}".SoUmiDedupRemSec.bam" | sort -k1,1 - > $tmp1 
grep -F -f <(sort -k1,1 ${readIDDir}${ID}".readIDs.txt" ) $tmp1 > $tmp2
cat $headerDir${ID}".SoUmiDedupRemSec.header" $tmp2 | samtools view -b | samtools sort > ${OUTBAMDIR}${ID}".snoRNAdownstream.read1.bam"

echo `wc -l ${readIDDir}${ID}.readIDs.txt`
echo `wc -l ${OUTBAMDIR}${ID}.snoRNAdownstream.read1.bam`

done