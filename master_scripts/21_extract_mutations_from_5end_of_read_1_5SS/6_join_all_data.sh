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

CIGAR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/3_5end_of_read1_bed_overlaps/"
REFNT="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/4_get_reference_nucleotide_of_5end_of_read_1_bed/"
SEQREAD="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/5_get_seqread_nucleotide_of_5end_of_read_1_bed/"
RELDIST="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/6_rel_dist_from_annotation/"

OUTDIR="11_extract_mutations_from_5end_of_read_1/1_5splice_sites/7_join_all_data/"

#make dirs if not already present
mkdir -p ${OUTDIR}


TMP1=${OUTDIR}"tmp1"
TMP2=${OUTDIR}"tmp2"


for INBED in ${RELDIST}*".rel_dist_to_anno.tab"
do 

#edit filename to save ID
ID=$(basename ${INBED} | cut -d "." -f 1)


join -1 4 -2 1 <(sed 's|/1\t|\t|g' ${INBED} | cut -f 1-6,10,13 | sort -k4,4) <(sed 's|/1\t|\t|g' ${CIGAR}${ID}".overlapping_annotation.5endofread1.bed" | cut -f 4,7,8 | sort -k1,1 ) |  tr " " "\t" > ${TMP1}
echo "first $ID"
join -1 1 -2 1 <(sort -k1,1 ${TMP1}) <(sort -k1,1 ${REFNT}${ID}".1stto3rdrefNt.tab") | tr " " "\t" > ${TMP2}
echo "second $ID"
join -1 1 -2 1 <(sort -k1,1 ${TMP2}) <(sort -k1,1 ${SEQREAD}${ID}".seq_read_nt_of_5end_of_read_1.tab") | tr " " "\t" > ${OUTDIR}${ID}".rel_dist_cigar_3refNT_3readNT_5SS.tab"
echo "third $ID"
rm ${TMP1} ${TMP2}

done

#remove final file if it exits 
rm -f ${OUTDIR}"xiCLIP_all.rel_dist_cigar_3refNT_3readNT_5SS.tab"

##cat all files together and append name to first feild

for INTAB in ${OUTDIR}*".rel_dist_cigar_3refNT_3readNT_5SS.tab"
do
ID=$(basename $INTAB | cut -d "." -f 1)

awk -v ID=$ID '{OFS="\t"}{print ID,$0}' ${INTAB} >> ${OUTDIR}"xiCLIP_all.rel_dist_cigar_3refNT_3readNT_5SS.tab"

done


