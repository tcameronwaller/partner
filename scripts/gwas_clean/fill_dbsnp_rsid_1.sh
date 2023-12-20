#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 December 2023
# Date, last execution: 19 December 2023
# Review: 19 December 2023
################################################################################
# Notes:



################################################################################
# Organize arguments.

path_file_reference_dbsnp=${1} #
path_directory_temporary=${2} #
path_bcftools=${3} #
report=${4} #

################################################################################
# Organize paths.

# Temporary files.
path_ftemp_dbsnp_biallelic_sites="${path_directory_temporary}/dbsnp_biallelic_sites.vcf.gz"
path_ftemp_dbsnp_extraction="${path_directory_temporary}/dbsnp_extraction.txt"
path_ftemp_dbsnp_extraction_alt_ref="${path_directory_temporary}/dbsnp_extraction_alt_ref.txt"
path_ftemp_dbsnp_extraction_ref_alt="${path_directory_temporary}/dbsnp_extraction_ref_alt.txt"


###########################################################################
# Execute procedure.

##########
# 1. Extract rsIDs and genomic coordinates from dbSNP.

# View header of dbSNP in Variant Call Format (VCF) to determine available tags
# for extraction.
#$path_bcftools head $path_file_reference_dbsnp
# "##INFO=<ID=RS,Number=1,Type=Integer,Description='dbSNP ID (i.e. rs number)'>"

# View records within 1000 Genomes file.
#$path_bcftools view --no-header --regions 7:70000-70010 $path_file_reference_dbsnp
# 7     70000 rs1041172676     G     A,C     ...     RS=1041172676;dbSNPBuildID=155;...
# 7     70010 rs1282455439     A     C       ...     RS=1282455439;dbSNPBuildID=155;...

# For simplicity, filter to sites (loci) that only have two allelic variants
# (biallelic sites).
#$path_bcftools norm --multiallelics +snps $path_file_reference_1kg_vcf | $path_bcftools view --no-header --min-alleles 2 --max-alleles 2 --types snps | head -10
$path_bcftools norm --multiallelics +snps $path_file_reference_dbsnp | $path_bcftools view --min-alleles 2 --max-alleles 2 --types snps > $path_ftemp_dbsnp_biallelic_sites

# Extract relevant information to a flat text table.
# https://samtools.github.io/bcftools/howtos/query.html
#$path_bcftools query -f '%ID %CHROM %POS %REF %ALT %INFO/EUR\n' $path_ftemp_dbsnp_biallelic_sites | head -10
echo "ID CHROM POS ALT REF RS_ID" > $path_ftemp_dbsnp_extraction
$path_bcftools query -f '%ID %CHROM %POS %ALT %REF %INFO/RS\n' $path_ftemp_dbsnp_biallelic_sites >> $path_ftemp_dbsnp_extraction



##########
# 2.1. For dbSNP, assemble information with unique identifiers specific to site
# (chromosome, position) and both reference and alternate alleles.

# Identifier format: <CHROM>_<POS>_<ALT>_<REF>
echo "identifier_merge ID CHROM POS ALT REF RS_ID" > $path_ftemp_dbsnp_extraction_alt_ref
cat $path_ftemp_dbsnp_extraction | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$4"_"$5), $1, $2, $3, $4, $5, $6
}' >> $path_ftemp_dbsnp_extraction_alt_ref

# Identifier format: <CHROM>_<POS>_<REF>_<ALT>
echo "identifier_merge ID CHROM POS ALT REF RS_ID" > $path_ftemp_dbsnp_extraction_ref_alt
cat $path_ftemp_dbsnp_extraction | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$5"_"$4), $1, $2, $3, $4, $5, $6
}' >> $path_ftemp_dbsnp_extraction_ref_alt



#
