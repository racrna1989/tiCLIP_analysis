#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

##################
#
#Purpose is to split bam files in to splice ($6 contains an N). and not split #Purpose - To extract read1 and ultimately extract the mutational profile of the reads 
#171220
#
##################

#activate environments 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#set wd
WD="/home/racrna/4_xiCLIP/"
cd ${WD}



MASTEROUTDIR="12_extract_splice_and_non_spliced_reads_close_to_exon_junctions/"
INDIR=${MASTEROUTDIR}"1_spliced_and_not_spliced_bams/"
OUTDIR=${MASTEROUTDIR}"2_counts_of_reads_overlapping_0nt_window_over_5end_or_3end_of_exons/"

ANNODILEREGIONS="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_exon_numbered.exonEnds0ntupdown.bed"
#ANNOFILEEXONS="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS.bed"

mkdir -p ${OUTDIR}


#for read 1, we can do intersect with -s option as read one is on cis strand. 
for INBAM in ${INDIR}/*read1*
do 

ID=$(basename ${INBAM})
echo ${ID}

#bedtools intersect -c -s -b ${INBAM} -a <(sed 's/^chr//g;s/|/_/g' $ANNOFILEEXONS | sort -k1,1 -k2,2n ) > ${OUTDIR}${ID}".Exon.counts"
bedtools intersect -c -s -b ${INBAM} -a <(sed 's/^chr//g;s/|/_/g' $ANNODILEREGIONS | sort -k1,1 -k2,2n )  > ${OUTDIR}${ID}".exonEnds0updown.counts"

done

#for read 2, we can do intersect with -s option as read one is on trans strand. 

for INBAM in ${INDIR}/*read2*
do 
ID=$(basename ${INBAM})


echo ${ID}

#bedtools intersect -c -S -b ${INBAM} -a <(sed 's/^chr//g;s/|/_/g' $ANNOFILEEXONS | sort -k1,1 -k2,2n ) > ${OUTDIR}${ID}.Exon.counts
bedtools intersect -c -S -b ${INBAM} -a <(sed 's/^chr//g;s/|/_/g' $ANNODILEREGIONS | sort -k1,1 -k2,2n )  > ${OUTDIR}${ID}.exonEnds0updown.counts

done




#for f in ${OUTDIR}*Exon*counts; do ID=$(basename $f | sed 's/.counts//g;s/.bam//g;s/\./_/g;s/-/_/g') ; awk -v ID=$ID '{OFS="\t"} $7 > 0 {print ID,$0}' $f >> ${MASTEROUTDIR}"xiCLIP_all_Exon_splice.count" ; done &
for f in ${OUTDIR}*exonEnds0updown.counts; do ID=$(basename $f | sed 's/.counts//g;s/.bam//g;s/\./_/g;s/-/_/g') ; awk -v ID=$ID '{OFS="\t"} $7 > 0 {print ID,$0}' $f >> ${MASTEROUTDIR}"xiCLIP_all_spliceSites0updown.count" ; done &

wait