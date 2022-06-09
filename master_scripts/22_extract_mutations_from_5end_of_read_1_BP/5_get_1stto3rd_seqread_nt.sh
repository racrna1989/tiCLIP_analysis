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

#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#filenames and dirs req

INDIR="11_extract_mutations_from_5end_of_read_1/2_branchpoints/2_read1_bam_overlaps/"
OUTDIR="11_extract_mutations_from_5end_of_read_1/2_branchpoints/5_get_seqread_nucleotide_of_5end_of_read_1_bed/"

#make dirs if not already present
mkdir -p ${OUTDIR}

for INBAM in ${INDIR}*".overlapping_annotation.read1.bam"
do 

#edit filename to save ID
ID=$(basename ${INBAM} | cut -d "." -f 1)

TMPFILE_1=${OUTDIR}"tmp1.tmp"
TMPFILE_2=${OUTDIR}"tmp2.tmp"
TMPFILE_3=${OUTDIR}"tmp3.tmp"
TMPFILE_4=${OUTDIR}"tmp4.tmp"
OUTFILE=${OUTDIR}${ID}".seq_read_nt_of_5end_of_read_1.tab"

#get seq read nt from:
#positive strand 

samtools view -f 96 ${INBAM} | cut -f 1 > ${TMPFILE_1}
samtools view -f 96 ${INBAM} | cut -f 10 | cut -c 1 - > ${TMPFILE_2}
samtools view -f 96 ${INBAM} | cut -f 10 | cut -c 2 - > ${TMPFILE_3}
samtools view -f 96 ${INBAM} | cut -f 10 | cut -c 3 - > ${TMPFILE_4}

paste ${TMPFILE_1} ${TMPFILE_2} ${TMPFILE_3} ${TMPFILE_4} > ${OUTFILE}

#negative strand 

samtools view -f 80 ${INBAM} | cut -f 1 > ${TMPFILE_1}
samtools view -f 80 ${INBAM} | cut -f 10 | rev | cut -c 1 - > ${TMPFILE_2}
samtools view -f 80 ${INBAM} | cut -f 10 | rev | cut -c 2 - > ${TMPFILE_3}
samtools view -f 80 ${INBAM} | cut -f 10 | rev | cut -c 3 - > ${TMPFILE_4}

paste ${TMPFILE_1} ${TMPFILE_2} ${TMPFILE_3} ${TMPFILE_4} >> ${OUTFILE}

sort -k1,1 ${OUTFILE} > tmp 

mv tmp ${OUTFILE}

rm ${TMPFILE_1} ${TMPFILE_2} ${TMPFILE_3} ${TMPFILE_4}

done





