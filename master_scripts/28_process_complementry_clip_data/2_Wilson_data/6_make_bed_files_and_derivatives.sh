#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 10
#SBATCH -A xiCLIP
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

# make bedfiles

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP


#make new working directory and cd to it
WD="/home/racrna/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/2_Wilson_EJC_iCLIP/"
cd $WD

INDIR="5_QC_mapped/"
OUTDIR="6_bed_files_and_derivatives/"

mkdir -p ${OUTDIR}

for INBAM in ${INDIR}*SoUmiDedupRemSec.bam;
do

#reconstruct name 

SAMPLE=$(basename $INBAM | sed 's/.SoUmiDedupRemSec.bam//g')


#make bed file plus and rev strand using read read1 and save to read1 directory intermediates are saved in tmp file
#bitflag 67 = read paired, mapped in properpair, first in pair i.e read1

echo -e "making beds from:\t" $SAMPLE
echo -e "making READ1 beds from:\t" $SAMPLE

samtools view -@ 8 -b -F 4 $INBAM | bedtools bamtobed -split -i - | sort -k1,1 -k2,2n > ${OUTDIR}${SAMPLE}".read.bed" &

echo -e "completed READ1 beds from:\t" $SAMPLE

##make 5prime pos bedfile from read bedfile (flag options are read mapped to forward strand)

echo -e "making 5primepos beds from:\t" $SAMPLE

samtools view -@ 8 -b -F 16 $INBAM | bedtools bamtobed | awk '{FS=OFS="\t"} {print $1,$2-1,$2,$4,$5,$6}' | awk -v OFS='\t' '{if ($2 < 0) print $1,$2+1,$3+1,$4,$5,$6; else print $0;}' > $SAMPLE".temp.5primepos.fwd.bed" &
samtools view -@ 8 -b -f 16 $INBAM | bedtools bamtobed | awk '{FS=OFS="\t"} {print $1,$3,$3+1,$4,$5,$6}' | awk -v OFS='\t' '{if ($2 < 0) print $1,$2+1,$3+1,$4,$5,$6; else print $0;}' > $SAMPLE".temp.5primepos.rev.bed" &

echo -e "completed 5primepos beds from:\t" $SAMPLE

wait

cat $SAMPLE".temp.5primepos.rev.bed" $SAMPLE".temp.5primepos.fwd.bed"| sort -k1,1 -k2,2n > ${OUTDIR}${SAMPLE}".5primepos.bed"

rm $SAMPLE".temp.5primepos.rev.bed" $SAMPLE".temp.5primepos.fwd.bed"


##extract 3endofRead2       (flag option is exclude unmapped reads)  	
echo -e "making 3endOfRead2 beds from:\t" $SAMPLE

samtools view -@ 8 -F 4 -b $INBAM | bedtools bamtobed | get3primeFlank.bash - 0 0 > ${OUTDIR}$SAMPLE".3endOfRead2.bed"
        	
echo -e "completed 3endOfRead2 beds from:\t" $SAMPLE
        	
echo -e "complete making beds from:\t" $SAMPLE

done

