#!/bin/bash
#script for uniting bg values with annotation file 200113
#all dir paths should contain / at the end. 
#InputVariables

bedAnnotation=$1
inBedDir=$2
outDir=$3

echo "Merging bedgraph counts with $bedAnnotation"
echo "Using bedgraphs in $inBedDir"
echo "Saving outputs in $outDir"

genome="/home/racrna/faststorage/GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.shortChrNames.cleaned.chr.sizes"

#construct tempDir file path from bedAnnotation filepath
tmpDir=$(dirname $bedAnnotation)
tmpDir=$tmpDir"/tmp/"

#constructing output suffix for filename. 
outFileName=$(basename $bedAnnotation | sed 's/bed//g')


#bedAnnotation="annotationFiles/HepG2_ENCODE.3ssSkippedIRandSpliced.191218.100ntupdown.binN201.bed"
#outDir="retainedIntronsHeLa/v4-sameTranscript-intronIDfixed-ENCODEspliceFactors/HepG2/200110_deconstructedSF3bandU2AF/"
#tmpDir="annotationFiles/tmp/"

#split annotation file and store paths as variables

mkdir -p $tmpDir

splitAnnOut=$(basename $bedAnnotation | sed 's/.bed//g')

if [ ! -s ${tmpDir}${splitAnnOut}"_plus.bed" ]; then
awk '{OFS="\t"} $6 =="+"' $bedAnnotation | sort -k1,1 -k2,2n > ${tmpDir}${splitAnnOut}"_plus.bed" &
awk '{OFS="\t"} $6 =="-"' $bedAnnotation | sort -k1,1 -k2,2n > ${tmpDir}${splitAnnOut}"_neg.bed" &
fi
wait

#save split annotation filepaths to appropriate variable
bedAnnotationPlus="${tmpDir}${splitAnnOut}_plus.bed"
bedAnnotationNeg="${tmpDir}${splitAnnOut}_neg.bed"

#these annotation files contain exons that will be substracted from the counts 
Exons="/home/racrna/faststorage/GRCh38/HeLa/custom/hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.bed"
snRNAs="/home/racrna/faststorage/GRCh38/snRNA/RNU_RNV_extentions.refFlat.bed"
snoRNAs="/home/racrna/faststorage/GRCh38/snoRNA/snoRNAs.GRCh38andrefGene.mature.bed"
#mk output dir

mkdir -p $outDir

#preform analysis

for inBedgraph in ${inBedDir}*rRNAScaled.fwd.bedgraph;
	do
	#make variables
	inBgPlus=$inBedgraph
	inBgNeg=$(echo $inBedgraph | sed 's/fwd.bedgraph/rev.bedgraph/')
	ID=$(basename $inBedgraph | sed 's/fwd.bedgraph//g' )
	RandPlus=$(($RANDOM))
	RandNeg=$(($RANDOM))

	echo "starting analysis for $ID but with the addage of removing exonic reads from intronic read density (eg small ncRNAs embedded in introns)"
	
	bedtools map -a $bedAnnotationPlus -b <( sort -k1,1 -k2,2n $inBgPlus | bedtools intersect -a - -b <(cat $Exons $snRNAs $snoRNAs | awk '{OFS="\t"} $6 =="+"' | sed 's/^chr//g' | sort -k1,1 -k2,2n ) ) -g $genome -c 4 -o sum -null 0 > /tmp/${ID}"_"${RandPlus}"_"${outFileName}"_tmp_plus" &
	
	
	bedtools map -a $bedAnnotationNeg -b <( sort -k1,1 -k2,2n $inBgNeg | bedtools intersect -a - -b <(cat $Exons $snRNAs $snoRNAs | awk '{OFS="\t"} $6 =="-"' | sed 's/^chr//g' | sort -k1,1 -k2,2n ) ) -g $genome -c 4 -o sum -null 0 > /tmp/${ID}"_"${RandNeg}"_"${outFileName}"_tmp_neg" &
		
	wait
	
	name=$(echo $ID | sed 's/\./\t/g' | cut -f 1 )
	
	cat /tmp/${ID}"_"${RandNeg}"_"${outFileName}"_tmp_neg" /tmp/${ID}"_"${RandPlus}"_"${outFileName}"_tmp_plus" | awk -v name=$name '{OFS="\t"}{print name, $0}' >> ${outDir}"all_Exons."${outFileName}"sense.counts"
	
	wait

	rm /tmp/${ID}"_"${RandNeg}"_"${outFileName}"_tmp_neg" /tmp/${ID}"_"${RandPlus}"_"${outFileName}"_tmp_plus"
	
	echo "completed coverage counts for $ID"
	
done


