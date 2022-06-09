#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP


WD="/home/racrna/4_xiCLIP/"
cd $WD

INTRONANNO="annotationFiles/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.bed"
INTRON_REF_FLAT_ANNO="/home/racrna/faststorage/GRCh38/refGene/introns.merged.exonssubstracted.annotated.bed"
SNORNAANNO="/home/racrna/faststorage/GRCh38/snoRNA/snoRNAs.GRCh38andrefGene.mature.bed"

OUTANNO="annotationFiles/introns_containing_snoRNAs_refFlat_and_Soren.bed"
SNORNA3END="annotationFiles/snoRNAs.GRCh38andrefGene.mature.3end.100ntupdown.binN201.bed"
SNORNA5END="annotationFiles/snoRNAs.GRCh38andrefGene.mature.5end.100ntupdown.binN201.bed"



# this part is used to fish out the introns that overlap with mature snoRNAs first using sorens annotation,
bedtools intersect -loj -s -a $SNORNAANNO -b <(sed 's/^chr//g' $INTRONANNO | sort ) | awk '{OFS="\t"}$8>-1{print $0}' > /tmp/snoRNA1.tmp
#and then using the refflat annotation, the final tmp files only contain those introns that overlap with snoRNAs, and it has them left joined.  
bedtools intersect -v -s -a $SNORNAANNO -b /tmp/snoRNA1.tmp |  bedtools intersect -loj -s -a - -b <(sed 's/^chr//g' $INTRON_REF_FLAT_ANNO | sort -k1,1 -k2,2n ) | awk '{OFS="\t"}$8>-1{print $0}' > /tmp/snoRNA2.tmp 

#finally, both these files are cat together and columns rearrange so geneID columns are spliced together
cat /tmp/snoRNA1.tmp /tmp/snoRNA2.tmp | awk '{OFS="\t"}{print $7,$8,$9,$4":::"$10,".",$12}' > $OUTANNO

OUTFILE=$(echo $OUTANNO | sed 's/.bed//g;s/\./_/g')
get5primeFlank.bash $OUTANNO 100 100 > $OUTFILE".100nt_5prime.bed"
get3primeFlank.bash $OUTANNO 100 100 > $OUTFILE".100nt_3prime.bed"
getWindows_3end.v3.sh $OUTANNO 100 100
getWindows_5end.v3.sh $OUTANNO 100 100

mv *binN201* annotationFiles/

#this labels the introns ends and then cats in the mature snoRNA file ends. 
cat <(awk '{OFS="\t"}{print $1,$2,$3,$4":::intron3ss",$5,$6}' $OUTFILE".3end.100ntupdown.binN201.bed") \
<(awk '{OFS="\t"}{print $1,$2,$3,$4":::intron5ss",$5,$6}' $OUTFILE".5end.100ntupdown.binN201.bed") \
$SNORNA3END \
$SNORNA5END \
| sort -k1,1 -k2,2n | grep -v random > "annotationFiles/hg38_snoRNAs_and_host_introns.100ntupdown.binN201.bed"

rm /tmp/snoRNA1.tmp /tmp/snoRNA2.tmp

#create index file for host gene, this can be used to normalise gene expression with sorens vector 
bedtools intersect -s -loj -a $SNORNAANNO -b <(sed 's/^chr//g' hg38_HeLa_trimmed_loci_major_primary_transcript.TotalExons.bed ) | \
awk '{OFS="\t"}$10 != "."{print $4,$10}' > snoRNA_containing_transcripts_index.tab

-a snoRNAs.GRCh38andrefGene.mature.bed | bedtools intersect -v -s -a - -b ~/faststorage/GRCh38/refGene/introns.merged.exonssubstracted.annotated.bed | grep -v nonintronic