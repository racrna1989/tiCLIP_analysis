#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 30G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xiCLIP


WD="/home/racrna/4_xiCLIP/"
cd $WD


INDIR="11_extract_mutations_from_5end_of_read_1/3_all_read_1_mut_comparisons/1_read_1_bams/"
GENOME="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.chr.sizes"



rRNAFACTOR="4_geneCount/rRNAFactor.tab"


for INBAM in ${INDIR}*mutation_at_first_nt.bam; do 

#this is to wrangle sample name to enable selection of scaling factor from rRNAFACTOR list
ID=$(basename $INBAM | cut -d "." -f 1)
	
SCALE=$(grep $ID $rRNAFACTOR | cut -f 2) 


#make bedgraph file NOTE grep -v K and G removes random chromosomes 
bedtools genomecov -5 -bg -scale ${SCALE} -strand + -ibam $INBAM  | grep -v "G" | grep -v "K" > ${INBAM/.bam/.rRNAScaled.plus.bedgraph} &
bedtools genomecov -5 -bg -scale ${SCALE} -strand - -ibam $INBAM  | grep -v "G" | grep -v "K" > ${INBAM/.bam/.rRNAScaled.minus.bedgraph} &

wait

bedGraphToBigWig ${INBAM/.bam/.rRNAScaled.plus.bedgraph} $GENOME ${INBAM/.bam/.rRNAScaled.plus.bw} &
bedGraphToBigWig ${INBAM/.bam/.rRNAScaled.minus.bedgraph} $GENOME ${INBAM/.bam/.rRNAScaled.minus.bw} &
wait


done