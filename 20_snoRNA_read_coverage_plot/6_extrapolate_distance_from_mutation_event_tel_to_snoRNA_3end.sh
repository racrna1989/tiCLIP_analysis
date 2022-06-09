#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 8
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

OUTDIR="10_snoRNA_read_coverage/4_getMutProfile/"
mkdir -p $OUTDIR

#filenames and dirs req


ANNODOWNSTREAMSNO="annotationFiles/snoRNAs.GRCh38andrefGene.maturetodowntreamexon.bed"
ANNOMATURESNO="annotationFiles/snoRNAs.GRCh38andrefGene.mature.bed"
regionBed="annotationFiles/snoRNAs.GRCh38andrefGene.maturetodowntreamexon.bed"
genomeFasta="~/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.fa"
BAMDIR="3_mapping/2_QC_mapped/"
headerDir=${OUTDIR}"2_samHeaders/"
OUTBAMDIR=${OUTDIR}"3_read1-snoRNA-bams/"
readIDDir=${OUTDIR}"1_readIDList/"
INDIR="10_snoRNA_read_coverage/3_map_snoRNA_to_read/"

#########
#PART 3##
#########

#extrapolate mutation position on read from sequencing.



mutationDir=${OUTDIR}"4_read_CIGARinfo/"
outAnalysisDir=${OUTDIR}"5_mutations_rel_to_snoRNA3end/"

mkdir -p $outAnalysisDir


for INBED in ${INDIR}*.bed
do
	#edit filename to save ID
	ID=$(basename $INBED | cut -d "." -f 1)
	
	#extrapolate mutation event rel to snoRNA3end ref pos $3 of mut minus ref pos of start of read $7 plus rel pos of start of read to 3end of sno (for +ve strand )
	#also only print those reads that have a mutation or deletion
	join -1 3 -2 4 \
	<(awk 'NR > 1' $mutationDir$ID".CIGARinfo.tab" | sort -k3,3 ) \
	<(sort -k4,4 $INBED) | sed 's/ /\t/g' | \
	cat <(echo -e "read.ID\tchr.mut\tmutation.refPos\tmutation.readPos\tOp.Ref.Base\tread.mapQ\tread.chr\tread.start\tread.end\tmapq\tread.strand\tstart.relTosno3End\tstop.relTosno3End\tsnoRNAID") - >  $outAnalysisDir$ID".extended.DistTosnoRNA3end.all.tab" 
	
	
	
	#awk '{OFS="\t"}{
	#if($10 == "+")
	#{print $0, (($3-$7)+$11)}
	#else if ($10 == "-")
	#{print $0, (($8-$3)+$12)}
	#}' | awk '{OFS="\t"}($5 !~ "N"){print $0}' > $outAnalysisDir$ID".extended.DistTosnoRNA3end.all.tab"
	
#	awk '{OFS="\t"}($5 !~ "S" && $5 !~ "N"){print $0}' $outAnalysisDir$ID".extended.DistTosnoRNA3end.all.tab" > $outAnalysisDir$ID".extended.DistTosnoRNA3end.MutandDel.tab"
	

done
