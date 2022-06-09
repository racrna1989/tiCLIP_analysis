#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 10
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00



source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#make new working directory and cd to it
WD="/home/racrna/4_xiCLIP/3_mapping"
mkdir -p $WD
cd $WD

#variables that contain filepaths for input fq, output bam and output qc'ed bam
FQINDIR="/home/racrna/4_xiCLIP/2_demultiplexing/4_RT_L3_demultiplexed_trimGalore/"
OUTBAMDIR="2_mapped/"
OUTMAPQCDIR="3_QC_mapped/"
mkdir -p $OUTBAMDIR $OUTMAPQCDIR


#mapping with no_softClipping
reference=/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa
genome=/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned
GTF=/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.77.cleaned.sorted.2.gtf
ss=/home/racrna/faststorage/GRCh38/genome/test/Homo_sapiens.GRCh38.77.cleaned.sorted.2.ss

#for loop to identify read 1 fq (BC included to exclude other file)
for FASTQ_R1_FILEPATH in $FQINDIR*CBP20*R1*BC*val*fq
do
##reconstruct flepath for read 2 & make name for out Sam
FASTQ_R2_FILEPATH=$(echo $FASTQ_R1_FILEPATH | sed 's/R1/R2/g;s/1_val_1/2_val_2/g')
OUTSAM=$(basename $FASTQ_R1_FILEPATH | sed 's/\./\t/g;s/_/\t/g;s/-/\t/g' | cut -f 2,3,10  | sed 's/\t/_/g')

echo -e "input filepath for read 1:\t"$FASTQ_R1_FILEPATH
echo -e "input filepath for read 2:\t"$FASTQ_R2_FILEPATH
echo -e "sample name \t" $OUTSAM
echo -e "beginning mapping of $OUTSAM"

##hisat2 new parameters, no-softclip. 
hisat2 \
-t \
--rna-strandness FR \
--fr \
--phred33 \
-p 10 \
--new-summary \
--no-softclip \
--known-splicesite-infile $ss \
-x $genome \
-1 $FASTQ_R1_FILEPATH \
-2 $FASTQ_R2_FILEPATH \
-S $OUTBAMDIR$OUTSAM".sam"

echo "beginning QC of mapped $OUTSAM"

##sort sam to bam
samtools import $reference $OUTBAMDIR$OUTSAM".sam" - | samtools sort -@ 10 - -o $OUTMAPQCDIR$OUTSAM".So.bam" 
##index bam for dedup
samtools index -@ 10 $OUTMAPQCDIR$OUTSAM".So.bam" 
##remove duplicates based on umi in header 
umi_tools dedup -I $OUTMAPQCDIR$OUTSAM".So.bam" --paired --output-stats=$OUTMAPQCDIR$OUTSAM".dedupStats." -S $OUTMAPQCDIR$OUTSAM".SoUmiDedup.bam" ;
##remove secondary alignments
samtools view -b -F256 $OUTMAPQCDIR$OUTSAM".SoUmiDedup.bam" > $OUTMAPQCDIR$OUTSAM".SoUmiDedupRemSec.bam"

done

