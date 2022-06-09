#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 10
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP

CHR_SIZES="/home/racrna/xiCLIP/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.chr.sizes"


#make new working directory and cd to it
WD="/home/racrna/4_xiCLIP/14_process_other_ALYREF_and_EJC_CLIP_data/2_Wilson_EJC_iCLIP/"
cd $WD

##start by counting number of reads in bed files

INDIR="5_QC_mapped/"
OUTDIR="7_mapped_reads/"

mkdir -p ${OUTDIR}

for INBAM in ${INDIR}*"SoUmiDedupRemSec.bam" ;  
	do 
	
	ID=$(basename $INBAM | cut -d"." -f 1)
	COUNT=$(samtools view -c $INBAM) 
	echo -e $ID"\t"$COUNT >> ${OUTDIR}"mapped_reads.tab"

done

#make bedgraphs and big wigs normalised to RPM

INDIR="6_bed_files_and_derivatives/"
OUTDIR="8_normalised_bigwigs_and_bedgraphs/"
COUNTFILE="7_mapped_reads/mapped_reads.tab"

mkdir -p ${OUTDIR}

for BED_FILEPATH in ${INDIR}*".bed"; 
	do
	#grep ID from list and compute scaling factor from reads mapped. 
	ID=$(basename $BED_FILEPATH | cut -d "." -f 1)
	COUNT=$(grep ${ID} ${COUNTFILE} | cut -f 2)
	RPM="1000000"
	SCALE=$(echo "scale=20; $RPM / $COUNT" | bc )


	bedtools genomecov -bg -split -scale $SCALE -strand + -i $BED_FILEPATH -g $CHR_SIZES | grep -v "G" | grep -v "K" > ${OUTDIR}${ID}".RPM.fwd.bedgraph"
	bedtools genomecov -bg -split -scale $SCALE -strand - -i $BED_FILEPATH -g $CHR_SIZES | grep -v "G" | grep -v "K" > ${OUTDIR}${ID}".RPM.rev.bedgraph"
	wait
	
	bedGraphToBigWig ${OUTDIR}${ID}".RPM.fwd.bedgraph" $CHR_SIZES ${OUTDIR}${ID}".RPM.fwd.bw" 
	bedGraphToBigWig ${OUTDIR}${ID}".RPM.rev.bedgraph" $CHR_SIZES ${OUTDIR}${ID}".RPM.rev.bw"
	wait
	
done