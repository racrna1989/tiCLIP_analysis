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

INTRON="/home/racrna/faststorage/GRCh38/HeLa/custom/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.bed"
ANNOBED="/home/racrna/faststorage/GRCh38/HeLa/custom/hg38_HeLa_annotated_mRNA_intronExon.bed"
OUTDIR1="4_geneCount/intron_exon_counts/"
OUTDIR2="4_geneCount/premRNA_vs_mRNA/"

INDIR="5_bedGraphs_and_derivatives/read2/"

mkdir -p $OUTDIR1 $OUTDIR2

for INBED in ${INDIR}*read2.bed 
do
	OUTNAME=$(basename ${INBED} | sed 's/.bed//g')
	echo "processing :" ${OUTNAME}
	#count overlaps between introns and exons, total overlaps may be more than total reads in bed
	bedtools intersect -c -s -a ${ANNOBED} -b ${INBED} > ${OUTDIR1}${OUTNAME}"Exon_Intron.counts"
	
	#exonic counts (read must reside within Exon)
	#INBED reads that overlap INTRON file are excluded, then counted with their overlap against exons 
	#Excluding reads that only map to Exons 
	
	bedtools intersect -c -s -a <(grep ":::intron" ${ANNOBED} ) -b <(bedtools intersect -s -a ${INBED} -b <(sed 's/^chr//g' ${INTRON} | sort -k1,1 -k2,2n)) |  awk -v name=$OUTNAME '{OFS="\t"}{print name"_premRNA", $0}' > ${OUTDIR2}${OUTNAME}"_premRNA.counts"
	
	bedtools intersect -c -s -a <(grep -v ":::intron" ${ANNOBED}) -b <(bedtools intersect -s -v -a ${INBED} -b <(sed 's/^chr//g' ${INTRON} | sort -k1,1 -k2,2n)) |  awk -v name=$OUTNAME '{OFS="\t"}{print name"_matureMRNA", $0}' > ${OUTDIR2}${OUTNAME}"_mature_mRNA.counts"

	#count total reads
	wc -l ${INBED} | cut -d " " -f 1 | awk -v name=$OUTNAME '{OFS="\t"}{print name, $0}' >> ${OUTDIR1}"xiCLIP.read2.totalcounts.200402.tab" 
	wc -l ${INBED} | cut -d " " -f 1 | awk -v name=$OUTNAME '{OFS="\t"}{print name, $0}' >> ${OUTDIR2}"xiCLIP.read2.totalcounts.200402.tab" 

done

#aggregate files 
for INCOUNT in ${OUTDIR1}*counts
do
	OUTNAME=$(basename ${INCOUNT} | sed 's/.counts//g;s/\./_/g')
	awk -v name=$OUTNAME '{OFS="\t"}$7 > 0 {print name,$0}' $INCOUNT >> ${OUTDIR1}"xiCLIP_intronExon.200402.counts"
done

