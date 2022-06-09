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

INDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/2_dist_from_annotation/"
OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_getMutProfile/"

mkdir -p $OUTDIR

#filenames and dirs req


GENOMEFASTA="~/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa"
ANNOBED="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.5SS.bed"
BAMDIR="3_mapping/2_QC_mapped/"

readIDDir=${OUTDIR}"1_readIDList/"
headerDir=${OUTDIR}"2_samHeaders/"
OUTBAMDIR=${OUTDIR}"3_read1_bam_overlaps/"

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
awk '{print $4}' ${INBED} | sort -k1,1 | sed 's|/1||g' | uniq > ${readIDDir}${ID}".readIDs.txt"
echo "printing header for read 1 $ID"	
#save header from original read1 bamfile.
samtools view -f 64 -H ${BAMDIR}${ID}".SoUmiDedupRemSec.bam" > $headerDir${ID}".SoUmiDedupRemSec.header" 


#mktmp file
TMPFILE_1=${OUTDIR}"tmp_1.tmp"
TMPFILE_2=${OUTDIR}"tmp_2.tmp"

echo "grepping IDs from sam"
#extract read1 of reads that overlap bam region file. Bame file region is annotation extended by 10nt either side. Grep reads, cat header and extracted reads
samtools view -f 64 -L <( awk '{OFS="\t"}{$2=($2-10); $3=($3+10); print $0}' ${ANNOBED}) ${BAMDIR}${ID}".SoUmiDedupRemSec.bam" | sort -k1,1 - > $TMPFILE_1 

grep -F -f <(sort -k1,1 ${readIDDir}${ID}".readIDs.txt" ) $TMPFILE_1 > $TMPFILE_2

cat $headerDir${ID}".SoUmiDedupRemSec.header" $TMPFILE_2 | samtools view -b | samtools sort - > ${OUTBAMDIR}${ID}".overlapping_annotation.read1.bam"

#note this count wont equal the readIDs in vs reads out as the samtools view -L option does not account for strand specificity. 
echo `wc -l ${readIDDir}${ID}.readIDs.txt`
echo `wc -l ${TMPFILE_2}`

rm $TMPFILE_1 $TMPFILE_2
done

