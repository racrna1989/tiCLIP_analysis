#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 20G
#SBATCH -c 1
#SBATCH -A xiCLIP
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ross.cordiner@mbg.au.dk

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rstudio

WD="/home/racrna/4_xiCLIP/"

cd $WD

inDir="4_geneCount/binnedGeneCounts-10kb/"

outDir="7_PCA_plots/"

mkdir -p $outDir


if [ -s ${outDir}"all.hg38HeLaSoren10kbins.NormTorRNAFactor.counts" ]; 
	then 
  		echo " ${outDir}all.hg38HeLaSoren10kbins.NormTorRNAFactor.counts exists and is not empty "
	else  

		for inCount in ${inDir}*_hg38HeLaSoren10kbins.NormTorRNAFactor.counts ; do
		
		cat $inCount >> ${outDir}"all.hg38HeLaSoren10kbins.NormTorRNAFactor.counts"

		done

fi

#run Rscript to process the files into a more managa

Rscript /home/racrna/4_xiCLIP/master_scripts/8_make_PCA/wrangle_data_file.R ${outDir}"all.hg38HeLaSoren10kbins.NormTorRNAFactor.counts" ${outDir}"all.hg38HeLaSoren10kbins.NormTorRNAFactor.tab"

