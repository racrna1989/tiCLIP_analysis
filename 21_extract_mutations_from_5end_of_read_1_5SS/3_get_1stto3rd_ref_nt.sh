#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH -c 1
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

#filenames and dirs req

GENOME="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa" 
INDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_5end_of_read1_bed_overlaps/"
OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/4_get_reference_nucleotide_of_5end_of_read_1_bed/"
GENOMESIZE="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.chr.sizes"

#make dirs if not already present
mkdir -p ${OUTDIR}

for INBED in ${INDIR}*".overlapping_annotation.5endofread1.bed"; 
do 

#edit filename to save ID
ID=$(basename ${INBED} | cut -d "." -f 1)

#get reference nucleotide
bedtools getfasta -s -tab -name -bed <(cut -f 1-6 ${INBED}) -fi ${GENOME} | sed 's|/1(.)||g' - > ${OUTDIR}${ID}".1strefNt.tab"

bedtools shift -i <(cut -f 1-6 ${INBED}) -g ${GENOMESIZE}  -p 1 -m -1 | bedtools getfasta -s -tab -name -bed - -fi ${GENOME} | sed 's|/1(.)\t|\t|g' - > ${OUTDIR}${ID}".2ndrefNt.tab"

bedtools shift -i <(cut -f 1-6 ${INBED}) -g ${GENOMESIZE}  -p 2 -m -2 | bedtools getfasta -s -tab -name -bed - -fi ${GENOME} | sed 's|/1(.)\t|\t|g' - > ${OUTDIR}${ID}".3rdrefNt.tab"


paste ${OUTDIR}${ID}".1strefNt.tab" <(cut -f 2 ${OUTDIR}${ID}".2ndrefNt.tab" ) <(cut -f 2 ${OUTDIR}${ID}".3rdrefNt.tab" )> ${OUTDIR}${ID}".1stto3rdrefNt.tab"

rm ${OUTDIR}${ID}".1strefNt.tab" ${OUTDIR}${ID}".2ndrefNt.tab" ${OUTDIR}${ID}".3rdrefNt.tab"

done





