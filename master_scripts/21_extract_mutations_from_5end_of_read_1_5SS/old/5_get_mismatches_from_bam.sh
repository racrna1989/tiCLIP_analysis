#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 8
#SBATCH -A xiCLIP
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk


WD="/home/racrna/4_xiCLIP/"

cd $WD

INDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_getMutProfile/3_read1_bam_overlaps/"
OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_getMutProfile/4_read_mismatches/"
mkdir -p $OUTDIR

#filenames and dirs req


GENOMEFASTA="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa"


#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

for INBAM in ${INDIR}*.bam
do
OUTNAME=$(basename $INBAM | sed 's/.bam//g')

samtools calmd -@ 8 -e $INBAM $GENOMEFASTA > ${OUTDIR}${OUTNAME}"_mismatch.sam"

grep -v "@" ${OUTDIR}${OUTNAME}"_mismatch.sam" > ${OUTDIR}${OUTNAME}"_mismatch.tab"


done


# samtools calmd -e 3_read1_bam_overlaps/RBM7_3_DMSO.overlapping_annotation.read1.bam $GENOMEFASTA | grep -v "@" | awk '{OFS="\t"}{print $1,$2,$5,$6,$9,$10,$(NF-4)}' | cut -f 7 | sort | cut -c 1-7 | uniq -c | sort -k1,1nr