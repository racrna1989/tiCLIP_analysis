#hg38_HeLa_Soren.10kbins.200409.bed
bed file of all TUs cut into 10kb segments.
$1 = chr
$2 = start
$3 = stop
$4 = TU info (TU name:::biotype:::totalExons:::TUsize)
$5 = segment number
$6 = strand

#hg38_HeLa_Soren.1kbins.200409.bed.zip
Same as above, except 1kb bins

#hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_exon_numbered.bed
bed file of all TUs exons
$1 = chr
$2 = start
$3 = stop
$4 = TU info (TU name:::biotype)
$5 = exon number in TU
$6 = strand

#hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_intron_numbered.bed
as above, except for introns


#hg38_HeLa_trimmed_loci_major_primary_isoform_annotated.exonNumber.sizeRange.TotalExonNumber.DistFromTSS.relDistToTSS_updatedWithUTRs_andSS.200429.bed
bed file containing same exons as "hg38_HeLa_trimmed_loci_major_primary_isoform_annotated_exon_numbered.bed", except with more information stored in $4
$1 = chr
$2 = start
$3 = stop
$4 = TU info (TU name:::biotype:::exonNumber:::totalExonsInTU:::exonSize:::cumSumExons*:::distToTSS**:::exonDesc***:::geneDesc****)
$5 = exon number in TU
$6 = strand

*cumSumExons = cumulativeDistanceToTSS. This is the sum of current and all previous exons.
**distToTSS =genomic distance to the TSS. This is the genomic distance of the 5' end of the exon to the 5' end of exon 1 in the TU. 
***majorExon = description of whether exon is predominately present in the transcript (can also be minor). 
****geneDesc = description of TU and exon position within TU. Can be singleExonicGene = monoExonic, multiExonicGene-firstExon, multiExonicGene-internalExon,multiExonicGene-lastExon
                               
#snoRNAs.GRCh38andrefGene.mature.bed
bed file of snoRNAs
$1 = chr
$2 = start
$3 = stop
$4 = snoRNA info (ENSG#:::layName:::biotype:::position[intronic/nonintronic])
$5 = .
$6 = strand
