#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

WD="/home/racrna/4_xiCLIP/"

cd $WD

INTRON="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.bed"
EXON="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_exon_numbered.bed"
OUTDIR="4_geneCount/intron_exon_counts/"
INDIR="5_bedGraphs_and_derivatives/read2/"
mkdir -p $OUTDIR


#process files and save in tmp directory 

ANNOBED="/tmp/exonIntron.bed"

cat $INTRON $EXON | sort -k1,1 -k2,2n | sed 's/^chr//g' > ${ANNOBED}
sed 's/^chr//g' $INTRON | sort -k1,1 -k2,2n > /tmp/intron.tmp
sed 's/^chr//g' $EXON | sort -k1,1 -k2,2n > /tmp/exon.tmp

EXON="/tmp/exon.tmp"
INTRON="/tmp/intron.tmp"


for INBED in ${INDIR}*read2.bed 
do
	OUTNAME=$(basename ${INBED} | sed 's/.bed//g;s/\./_/g')
	
	echo "processing :" ${OUTNAME}
	
	#count total reads
	wc -l ${INBED} | cut -d " " -f 1 | awk -v name=$OUTNAME '{OFS="\t"}{print name, $0}' >> ${OUTDIR}"xiCLIP.read2.totalcounts.200402.tab" &
	
	#count overlaps between introns and exons, total overlaps may be more than total reads in bed
	bedtools intersect -c -s -a ${ANNOBED} -b ${INBED} | \
	awk -v name=$OUTNAME '{OFS="\t"}$7 > 0 {print name,$0}' >> ${OUTDIR}"xiCLIP_intronExon.200820.counts" & 

	
	#count overlaps between introns and reads mapped exclusively to exons
	bedtools intersect -c -s -a ${INTRON} -b ${INBED} > ${OUTDIR}${OUTNAME}".Intron.counts" &
	bedtools intersect -c -s -a ${EXON} -b <(bedtools intersect -s -v -a ${INBED} -b ${INTRON}) > ${OUTDIR}${OUTNAME}".Exon.counts" & 
	
	wait
	cat ${OUTDIR}${OUTNAME}".Intron.counts" ${OUTDIR}${OUTNAME}".Exon.counts" | awk -v name=$OUTNAME '{OFS="\t"}$7 > 0 {print name,$0}' >> ${OUTDIR}"xiCLIP_intron_exclusiveExon.200820.counts"
	rm ${OUTDIR}${OUTNAME}".Intron.counts" ${OUTDIR}${OUTNAME}".Exon.counts"	
	
	

done

rm ${ANNOBED}
rm $(INTRON}
rm ${EXON}