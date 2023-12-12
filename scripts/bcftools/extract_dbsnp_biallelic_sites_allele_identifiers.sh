#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 12 December 2023
# Date, last execution: 12 December 2023
# Review: 12 December 2023
################################################################################
# Notes:

# This script extracts reference SNP cluster identifiers (rsIDs) for biallelic
# sites from the dbSNP reference and gives them chromosome-position-allele
# identifiers for merging to other tables.

# The script below is a companion to perform subsequent merging operations with
# the tables that this script creates.
# "/.../partner/scripts/gwas_clean/fill_reference_snp_cluster_identifier.sh"

################################################################################



################################################################################
# Organize paths.

# Tools.
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_reference_dbsnp="${path_directory_reference}/dbsnp/grch37_chromosome" # dbSNP build 155; accession: TCW; 2023-02-06; chromosome translation: TCW; 2023-02-15

# Files.
path_file_reference_dbsnp="${path_directory_reference_dbsnp}/GCF_000001405.25.gz" # dbSNP build 155; accession: TCW; 2023-02-06; chromosome translation: TCW; 2023-02-15
path_file_dbsnp_biallelic_sites="${path_directory_reference_dbsnp}/dbsnp_biallelic_sites.vcf.gz"
path_file_dbsnp_extraction="${path_directory_reference_dbsnp}/dbsnp_extraction.txt"
path_file_dbsnp_extraction_alt_ref="${path_directory_reference_dbsnp}/dbsnp_extraction_alt_ref.txt"
path_file_dbsnp_extraction_ref_alt="${path_directory_reference_dbsnp}/dbsnp_extraction_ref_alt.txt"

# Initialize files.
rm $path_file_dbsnp_biallelic_sites
rm $path_file_dbsnp_extraction
rm $path_file_dbsnp_extraction_alt_ref
rm $path_file_dbsnp_extraction_ref_alt

# Initialize directories.
cd $path_directory_reference_dbsnp



###########################################################################
# Parameters.

report="true"

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
$path_bcftools norm --multiallelics +snps $path_file_reference_dbsnp | $path_bcftools view --min-alleles 2 --max-alleles 2 --types snps > $path_file_dbsnp_biallelic_sites

# Extract relevant information to a flat text table.
# https://samtools.github.io/bcftools/howtos/query.html
#$path_bcftools query -f '%ID %CHROM %POS %REF %ALT %INFO/EUR\n' $path_file_dbsnp_biallelic_sites | head -10
echo "ID CHROM POS ALT REF RS_ID" > $path_file_dbsnp_extraction
$path_bcftools query -f '%ID %CHROM %POS %ALT %REF %INFO/RS\n' $path_file_dbsnp_biallelic_sites >> $path_file_dbsnp_extraction



##########
# 2. For dbSNP, assemble information with unique identifiers specific to site
# (chromosome, position) and both reference and alternate alleles.

# Identifier format: <CHROM>_<POS>_<ALT>_<REF>
echo "identifier_merge ID CHROM POS ALT REF RS_ID" > $path_file_dbsnp_extraction_alt_ref
cat $path_file_dbsnp_extraction | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$4"_"$5), $1, $2, $3, $4, $5, $6
}' >> $path_file_dbsnp_extraction_alt_ref

# Identifier format: <CHROM>_<POS>_<REF>_<ALT>
echo "identifier_merge ID CHROM POS ALT REF RS_ID" > $path_file_dbsnp_extraction_ref_alt
cat $path_file_dbsnp_extraction | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$5"_"$4), $1, $2, $3, $4, $5, $6
}' >> $path_file_dbsnp_extraction_ref_alt



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Extract information from dbSNP reference."
  echo "----------"
  echo "path to source dbSNP file: " $path_file_reference_dbsnp
  echo "path to product dbSNP extraction file: " $path_file_dbsnp_extraction
  echo "----------"
  echo "dbSNP extraction:"
  cat $path_file_dbsnp_extraction | head -5
  echo "- - Count of lines (rows, records) in table:"
  cat $path_file_dbsnp_extraction | wc -l
  echo "----------"
  echo "dbSNP extraction identifier: <CHROM>_<POS>_<ALT>_<REF>"
  cat $path_file_dbsnp_extraction_alt_ref | head -5
  echo "- - Count of lines (rows, records) in table:"
  cat $path_file_dbsnp_extraction_alt_ref | wc -l
  echo "----------"
  echo "dbSNP extraction identifier: <CHROM>_<POS>_<REF>_<ALT>"
  cat $path_file_dbsnp_extraction_ref_alt | head -5
  echo "- - Count of lines (rows, records) in table:"
  cat $path_file_dbsnp_extraction_ref_alt | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi



#
