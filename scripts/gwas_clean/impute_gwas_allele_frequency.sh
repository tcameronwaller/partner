#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 February 2023
# Date, last execution: 2 August 2023
# Review: TCW; 2 August 2023
################################################################################
# Note

# This script extracts from genomic variation (Variant Call Format; VCF) in
# Phase 3 of the 1000 Genomes Project frequencies of alternate alleles for the
# European ("EUR") superpopulation and introduces these to a set of GWAS summary
# statistics.
# For simplicity, this script filters to biallelic variant sites.

# GWAS summary statistics Format (Team Standard)
# Effect allele: "A1"
# Delimiter: white space
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

################################################################################




################################################################################
# Organize arguments.

path_file_gwas_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_gwas_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
report=${3} # whether to print reports

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_reference_1kg="${path_directory_reference}/1000_genomes_phase_3" # genomic variation (Variant Call Format; VCF); accession: TCW; ___
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="$(dirname $path_file_gwas_product)"
name_base_file_gwas_product="$(basename $path_file_gwas_product .txt.gz)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_gwas_product}" # must be unique

# Files.
path_file_reference_1kg_vcf="${path_directory_reference_1kg}/1000GENOMES-phase_3.vcf.gz" # accession: TCW; 2023-02-22;

# Temporary files.
path_ftemp_1kg_biallelic_sites="${path_directory_temporary}/1kg_biallelic_sites.vcf.gz"
path_ftemp_1kg_frequency="${path_directory_temporary}/1kg_european_allele_frequency.txt"
path_ftemp_1kg_frequency_identifier="${path_directory_temporary}/1kg_european_allele_frequency_identifier.txt"
path_ftemp_gwas_identifier="${path_directory_temporary}/gwas_source_identifier.txt"
path_ftemp_frequency_gwas_merge="${path_directory_temporary}/frequency_gwas_merge.txt"
path_ftemp_frequency_gwas_merge_filter="${path_directory_temporary}/frequency_gwas_merge_filter.txt"
path_ftemp_product_format="${path_directory_temporary}/gwas_product_format.txt"

# Initialize files.
rm $path_file_gwas_product

# Initialize directories.
#rm -r $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product

###########################################################################
# Execute procedure.

##########
# 1. Extract allelic frequencies for European ancestry from genomic variation in
#    Phase 3 of the 1000 Genomes Project.

# BCFTools
# A bioinformatic tool for tasks on genotype files in Variant Call Format (VCF).
# Documentation from SamTools (https://samtools.github.io/bcftools/howtos/index.html).
# Documentation manual (https://samtools.github.io/bcftools/bcftools.html).
# https://plink.readthedocs.io/en/latest/bcftools_mani/
# https://www.htslib.org/howtos/headers.html

# View header of 1000 Genomes file in Variant Call Format (VCF) to determine
# available tags for extraction.
#$path_bcftools head $path_file_reference_1kg_vcf
# "##INFO=<ID=EUR,Number=A,Type=Float,Description='Allele frequency for European populations in the 1000 Genomes Project Phase 3'>"

# View records within 1000 Genomes file.
#$path_bcftools view --no-header --regions 7:1000-1010 $path_file_reference_1kg_vcf

# Filter to sites (loci) that only have two allelic variants (biallelic sites).
#$path_bcftools norm --multiallelics +snps $path_file_reference_1kg_vcf | $path_bcftools view --no-header --min-alleles 2 --max-alleles 2 --types snps | head -10
$path_bcftools norm --multiallelics +snps $path_file_reference_1kg_vcf | $path_bcftools view --min-alleles 2 --max-alleles 2 --types snps > $path_ftemp_1kg_biallelic_sites

# Extract relevant information to a flat text table.
# https://samtools.github.io/bcftools/howtos/query.html
#$path_bcftools query -f '%ID %CHROM %POS %REF %ALT %INFO/EUR\n' $path_ftemp_1kg_biallelic_sites | head -10
echo "ID CHROM POS ALT REF AF_EUR" > $path_ftemp_1kg_frequency
$path_bcftools query -f '%ID %CHROM %POS %ALT %REF %INFO/EUR\n' $path_ftemp_1kg_biallelic_sites >> $path_ftemp_1kg_frequency

##########
# 2. Assemble unique identifiers specific to site (chromosome, position) and
#    allele (alternate allele).

# Alternate allele frequencies from 1000 Genomes.
echo "identifier_merge ID CHROM POS ALT REF AF_EUR" > $path_ftemp_1kg_frequency_identifier
cat $path_ftemp_1kg_frequency | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$4), $1, $2, $3, $4, $5, $6
}' >> $path_ftemp_1kg_frequency_identifier

# GWAS summary statistics.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_gwas_identifier
zcat $path_file_gwas_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$4), $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14
}' >> $path_ftemp_gwas_identifier

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Alternate allele frequencies from 1000 Genomes European ancestry."
  echo "- - table before merge:"
  head -30 $path_ftemp_1kg_frequency_identifier
  echo "----------"
  echo "----------"
  echo "----------"
  echo "GWAS summary statistics."
  echo "- - table before merge:"
  head -30 $path_ftemp_gwas_identifier
  echo "----------"
  echo "----------"
  echo "----------"
fi

##########
# 3. Introduce allelic frequencies to GWAS summary statistics.
# Merge text tables by a common identifier in "awk".
# Bash command "join" might also be capable of this type of merge or join, but
# it might require the same identifiers in sort order.
# For computational efficiency, the first file ought to be the subset of the
# second file.
# Table 1: 15 total columns with merge identifier in column 1.
# Table 2: 7 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "identifier_merge ID CHROM POS ALT REF AF_EUR SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_frequency_gwas_merge
awk 'FNR==NR{a[$1]=$2FS$3FS$4FS$5FS$6FS$7FS$8FS$9FS$10FS$11FS$12FS$13FS$14FS$15; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, a[$1]}
' $path_ftemp_gwas_identifier $path_ftemp_1kg_frequency_identifier > $path_ftemp_frequency_gwas_merge

cat $path_ftemp_frequency_gwas_merge | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_ftemp_frequency_gwas_merge_filter
cat $path_ftemp_frequency_gwas_merge | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 21)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( $1 == "NA" )
    # Empty record indicates that second file that did not match first file.
    next
  else if (($3 != $9) || ($4 != $10) || (toupper($5) != toupper($11)))
    # Chromosome, position, or alternate allele do not match.
    next
  else
    # Keep record from merge of first file and second file.
    print $0
  }' >> $path_ftemp_frequency_gwas_merge_filter

echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_product_format
cat $path_ftemp_frequency_gwas_merge_filter | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print $8, $9, $10, toupper($11), toupper($12), $7, $14, $15, $16, $17, $18, $19, $20, $21
}' >> $path_ftemp_product_format

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "GWAS summary statistics with alternate allele frequencies."
  echo "- - table after merge:"
  head -30 $path_ftemp_product_format
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Compress file format.
gzip -cvf $path_ftemp_product_format > $path_file_gwas_product



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Procedure complete."
  echo "- - Count of lines in source GWAS summary statistics:"
  zcat $path_file_gwas_source | wc -l
  echo "- - Count of lines in product GWAS summary statistics:"
  zcat $path_file_gwas_product | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi



##########
# Remove temporary directories and files.
# Suppress this block for debugging.
if true; then
  rm -r $path_directory_temporary
fi



#
