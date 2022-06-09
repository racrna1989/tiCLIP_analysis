#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 20
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk


WD="/home/racrna/4_xiCLIP/"

cd $WD

INDIR="3_mapping/2_QC_mapped/"
OUTDIR="11_extract_mutations_from_5end_of_read_1/3_all_read_1_mut_comparisons/1_read_1_bams/"

OUTDIR_2="11_extract_mutations_from_5end_of_read_1/3_all_read_1_mut_comparisons/2_read_1_sam_split_by_containing_mutation/"
mkdir -p $OUTDIR

#filenames and dirs req


GENOMEFASTA="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa"


#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

for INBAM in ${INDIR}*RBM7*SoUmiDedupRemSec.bam
do

ID=$(basename ${INBAM} | cut -d "." -f 1)

#samtools view -f 64 -b ${INBAM} > ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.bam"
#samtools view -f 64 -H ${INBAM} > ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.header"

#samtools calmd -@ 20 -e ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.bam" ${GENOMEFASTA} > ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.sam"


#####
#rm -f ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.mutation_at_first_nt.sam"


#grep reads that have a mismatch on first nucleotide in reads 
#fwd strand 96 = mate on reverse, first in pair

#samtools view -f 96 ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.sam" | grep -e $'\t''[ATGC]=' >> ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.mutation_at_first_nt.sam"

#rev strand 80 = read on reverse, first in pair

#samtools view -f 80 ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.sam" | grep -e '=[ATGC]'$'\t' >> ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.mutation_at_first_nt.sam"

#echo ${ID} "mutation" 

#cat ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.header" ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.mutation_at_first_nt.sam" | samtools view -b | samtools sort - > ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.mutation_at_first_nt.bam"


#############

#rm -f ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.no_mutation_at_first_nt.sam"

#grep reads that dont have a mismatch on first nucleotide in reads 
#fwd strand 96 = mate on reverse, first in pair

#samtools view -f 96 ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.sam" | grep -v -e $'\t''[ATGC]=' >> ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.no_mutation_at_first_nt.sam" 

#rev strand 80 = read on reverse, first in pair

#samtools view -f 80 ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.sam" | grep -v -e '=[ATGC]'$'\t' >> ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.no_mutation_at_first_nt.sam" 

#echo ${ID} "no mutation" 

#cat ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.header" ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.no_mutation_at_first_nt.sam" | samtools view -b | samtools sort - > ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.no_mutation_at_first_nt.bam"

###########


samtools index ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.no_mutation_at_first_nt.bam"
samtools index ${OUTDIR}${ID}".SoUmiDedupRemSec.read1.calmd.mutation_at_first_nt.bam"


done


# samtools calmd -e 3_read1_bam_overlaps/RBM7_3_DMSO.overlapping_annotation.read1.bam $GENOMEFASTA | grep -v "@" | awk '{OFS="\t"}{print $1,$2,$5,$6,$9,$10,$(NF-4)}' | cut -f 7 | sort | cut -c 1-7 | uniq -c | sort -k1,1nr