#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 10
#SBATCH -A xiCLIP
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk



source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

#filepaths
genome=~/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.chr.sizes


#make new working directory and cd to it
WD="/home/racrna/4_xiCLIP/"

cd $WD

inDir="3_mapping/2_QC_mapped/"
outDir="/home/racrna/4_xiCLIP/5_bedGraphs_and_derivatives/"

outDir5PRIMEPOS="5_bedGraphs_and_derivatives/5primepos/"
outDirREAD1="5_bedGraphs_and_derivatives/read1/"
outDirREAD2="5_bedGraphs_and_derivatives/read2/"
outDir3END="5_bedGraphs_and_derivatives/3endOfRead2/"

mkdir -p $outDir 
mkdir -p $outDir5PRIMEPOS $outDirREAD2 $outDirREAD1 $outDir3END

#programmes

#selecting 88 or 83, selects reads that are mapped to the positive or negative strand.
#99 = read paired (0x1), read mapped in proper pair (0x2), mate reverse strand (0x20), first in pair (0x40)
#83 = read paired (0x1), read mapped in proper pair (0x2), read reverse strand (0x10), first in pair (0x40)

for inBam in $inDir*ALYREF*.SoUmiDedupRemSec.bam
do
#reconstruct name 

SAMPLE=$(basename $inBam | sed 's/.SoUmiDedupRemSec.bam//g')


#make bed file plus and rev strand using read read1 and save to read1 directory intermediates are saved in tmp file
#bitflag 67 = read paired, mapped in properpair, first in pair i.e read1

echo -e "making beds from:\t" $SAMPLE
echo -e "making READ1 beds from:\t" $SAMPLE
samtools view -@ 8 -b -f 67 $inBam | bedtools bamtobed -split -i - | sort -k1,1 -k2,2n > $outDirREAD1$SAMPLE".read1.bed" &
echo -e "completed READ1 beds from:\t" $SAMPLE

##make 5prime pos bedfile from read1 bedfile

echo -e "making 5primepos beds from:\t" $SAMPLE

samtools view -@ 8 -b -f 99 $inBam | bedtools bamtobed | awk '{FS=OFS="\t"} {print $1,$2-1,$2,$4,$5,$6}' | awk -v OFS='\t' '{if ($2 < 0) print $1,$2+1,$3+1,$4,$5,$6; else print $0;}' > /tmp/$SAMPLE".temp.5primepos.fwd.bed" &
samtools view -@ 8 -b -f 83 $inBam | bedtools bamtobed | awk '{FS=OFS="\t"} {print $1,$3,$3+1,$4,$5,$6}' | awk -v OFS='\t' '{if ($2 < 0) print $1,$2+1,$3+1,$4,$5,$6; else print $0;}' > /tmp/$SAMPLE".temp.5primepos.rev.bed" &
echo -e "completed 5primepos beds from:\t" $SAMPLE

wait

cat /tmp/$SAMPLE".temp.5primepos.rev.bed" /tmp/$SAMPLE".temp.5primepos.fwd.bed"| sort -k1,1 -k2,2n > $outDir5PRIMEPOS$SAMPLE".5primepos.bed"

##extract read2 from file (same flag as above, except second in pair activated x80) also swaps the strands with awk script. The read 2 is split, to identify coverage over spliced exons.
##this file cannot be used to generate 3endOfread2 as it creates two new bed reads, with nonidentical 3'ends. 

echo -e "making READ2 beds from:\t" $SAMPLE

samtools view -@ 8 -b -f 131 $inBam | bedtools bamtobed -split | \
awk '{OFS="\t"}{
        	if ($6 =="+"){$6="-"; print $0}
        	else if ($6 == "-"){$6="+"; print $0}
        	}' > $outDirREAD2$SAMPLE".read2.bed" 
        	
echo -e "completed READ2 beds from:\t" $SAMPLE

##extract 3endofRead2        	
echo -e "making 3endOfRead2 beds from:\t" $SAMPLE

samtools view -@ 8 -f 131 -b $inBam | bedtools bamtobed |
        awk '{OFS="\t"}{
        	if ($6 =="+"){$6="-"; print $0}
        	else if ($6 == "-"){$6="+"; print $0}
        	}' | \
        	get3primeFlank.bash - 0 0 > $outDir3END$SAMPLE".3endOfRead2.bed"
        	
echo -e "completed 3endOfRead2 beds from:\t" $SAMPLE
        	
echo -e "complete making beds from:\t" $SAMPLE

done

