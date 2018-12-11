#! /usr/bin/env bash

########################################################################
infile=$1 ###vcf file to annotate
outFile=$2 ###prefix for the output files

########Usage vcfAnnotate.sh <input_vcf_file> <output_prefix>##########

python2.7 scripts/filtervcf_v5.py $infile $outFile
	
python2.7 scripts/ResPred_ann_v9.py scripts/"refdict_v5" scripts/"query_snps.txt" scripts/"query_indels.txt" scripts/"coord_v3" $outFile"_filtered.vcf" $outFile scripts/"gene_list" ./
