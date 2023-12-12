#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 December 2023
# Date, last execution: 12 December 2023
# Review: 12 December 2023
################################################################################
# Notes:

# This script extracts reference SNP cluster identifiers (rsIDs) from the dbSNP
# reference and matches them to genomic coordinates in GWAS summary statistics
# to fill missing rsIDs in the GWAS summary statistics.

# Allele designations.
# 1. In the dbSNP VCF file, the reference allele ("REF") corresponds to the
# reference human genome, and the alternate allele(s) ("ALT") correspond to the
# alternative allele(s) that sometimes replaces the reference allele at that
# position in the human genome.
# 2. In GWAS summary statistics, either the dbSNP reference allele or the dbSNP
# alternate allele might correspond to the effect allele. Whichever allele is
# the effect allele corresponds to the directionality of the GWAS association
# effect parameter (beta).

# The script below is a companion to perform the first few operations of this
# script to save time in subsequent iterations.
# "/.../partner/scripts/bcftools/extract_dbsnp_biallelic_sites_allele_identifiers.sh"


################################################################################


################################################################################
# Organize arguments.

path_file_gwas_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_gwas_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
report=${3} # whether to print reports


################################################################################
# Organize paths.

# Tools.
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_reference_dbsnp="${path_directory_reference}/dbsnp/grch37_chromosome" # dbSNP build 155; accession: TCW; 2023-02-06; chromosome translation: TCW; 2023-02-15
path_directory_product="$(dirname $path_file_gwas_product)"
name_base_file_gwas_product="$(basename $path_file_gwas_product .txt.gz)"
path_directory_temporary="${path_directory_product}/temporary_2713956_${name_base_file_gwas_product}" # must be unique

# Files.
path_file_reference_dbsnp="${path_directory_reference_dbsnp}/GCF_000001405.25.gz" # dbSNP build 155; accession: TCW; 2023-02-06; chromosome translation: TCW; 2023-02-15
if true; then
  path_ftemp_dbsnp_extraction="${path_directory_reference_dbsnp}/dbsnp_extraction.txt"
  path_ftemp_dbsnp_extraction_alt_ref="${path_directory_reference_dbsnp}/dbsnp_extraction_alt_ref.txt"
  path_ftemp_dbsnp_extraction_ref_alt="${path_directory_reference_dbsnp}/dbsnp_extraction_ref_alt.txt"
fi

# Temporary files.
if false; then
  path_ftemp_dbsnp_biallelic_sites="${path_directory_temporary}/dbsnp_biallelic_sites.vcf.gz"
  path_ftemp_dbsnp_extraction="${path_directory_temporary}/dbsnp_extraction.txt"
  path_ftemp_dbsnp_extraction_alt_ref="${path_directory_temporary}/dbsnp_extraction_alt_ref.txt"
  path_ftemp_dbsnp_extraction_ref_alt="${path_directory_temporary}/dbsnp_extraction_ref_alt.txt"
fi
path_ftemp_gwas_identifier_a1_a2="${path_directory_temporary}/gwas_source_identifier_a1_a2.txt"
path_ftemp_merge_alt_ref="${path_directory_temporary}/merge_alt_ref.txt"
path_ftemp_merge_alt_ref_clean="${path_directory_temporary}/merge_alt_ref_clean.txt"
path_ftemp_merge_ref_alt="${path_directory_temporary}/merge_ref_alt.txt"
path_ftemp_merge_ref_alt_clean="${path_directory_temporary}/merge_ref_alt_clean.txt"
path_ftemp_merge_priority="${path_directory_temporary}/merge_priority.txt"
path_ftemp_product_format="${path_directory_temporary}/gwas_product_format.txt"

# Initialize files.
rm $path_file_gwas_product

# Initialize directories.
rm -r $path_directory_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product



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

if false; then
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
fi



##########
# 3. For GWAS summary statistics, assemble information with unique identifiers
# specific to site (chromosome, position) and both effect and other alleles.

# GWAS summary statistics.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_gwas_identifier_a1_a2
zcat $path_file_gwas_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print ($2"_"$3"_"$4"_"$5), $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14
}' >> $path_ftemp_gwas_identifier_a1_a2



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "dbSNP extraction identifier: <CHROM>_<POS>_<ALT>_<REF>"
  echo "- - table before merge:"
  head -5 $path_ftemp_dbsnp_extraction_alt_ref
  echo "----------"
  echo "dbSNP extraction identifier: <CHROM>_<POS>_<REF>_<ALT>"
  echo "- - table before merge:"
  head -5 $path_ftemp_dbsnp_extraction_ref_alt
  echo "----------"
  echo "GWAS summary statistics identifier: <CHROM>_<POS>_<A1>_<A2>"
  echo "- - table before merge:"
  head -5 $path_ftemp_gwas_identifier_a1_a2
  echo "----------"
  echo "----------"
  echo "----------"
fi



##########
# 4. Merge dbSNP information to GWAS summary statistics.
# Merge text tables by a common identifier in "awk".
# Bash command "join" might also be capable of this type of merge or join, but
# it might require the same identifiers in sort order.
# For computational efficiency, the first file ought to be the subset of the
# second file.

# Merge 1.
# Table 1: GWAS summary statistics
# Table 1: 15 total columns with merge identifier in column 1.
# Table 2: dbSNP extraction with identifier format <CHROM>_<POS>_<ALT>_<REF>.
# Table 2: 7 total columns with merge identifier in column 1.
# Merged table: 21 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "identifier_merge ID CHROM POS ALT REF RS_ID SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_merge_alt_ref
awk 'FNR==NR{a[$1]=$2FS$3FS$4FS$5FS$6FS$7FS$8FS$9FS$10FS$11FS$12FS$13FS$14FS$15; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, a[$1]}
' $path_ftemp_gwas_identifier_a1_a2 $path_ftemp_dbsnp_extraction_alt_ref > $path_ftemp_merge_alt_ref

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 1:"
  head -10 $path_ftemp_merge_alt_ref
  echo "----------"
  echo "----------"
  echo "----------"
fi

# $1               $2 $3    $4  $5  $6  $7    $8  $9  $10 $11 $12 $13  $14  $15 $16 $17 $18 $19  $20   $21
# identifier_merge ID CHROM POS ALT REF RS_ID SNP CHR BP  A1  A2  A1AF BETA SE  P   N   Z   INFO NCASE NCONT

# Check and filter the information from the merge.
echo "identifier_merge ID_ar ALT_ar REF_ar SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_merge_alt_ref_clean
#cat $path_ftemp_merge_alt_ref | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_ftemp_merge_alt_ref_clean
cat $path_ftemp_merge_alt_ref | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 21)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( $1 == "NA" )
    # Missing identifier indicates that second file did not match first file.
    # Skip the record.
    next
  else if (($9 != $3) || ($10 != $4) || ((toupper($11) != toupper($5)) && (toupper($11) != toupper($6))))
    # Chromosome, position, or alleles do not match.
    # Notice the comparison of A1 to both ALT and REF.
    next
  else
    # Keep record from successful merge.
    print $1, $2, $5, $6, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21
}' >> $path_ftemp_merge_alt_ref_clean



# Merge 2.
# Table 1: Table from Merge 1.
# Table 1: 18 total columns with merge identifier in column 1.
# Table 2: dbSNP extraction with identifier format <CHROM>_<POS>_<REF>_<ALT>.
# Table 2: 7 total columns with merge identifier in column 1.
# Merged table: 24 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "identifier_merge ID CHROM POS ALT REF RS_ID ID_ar SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_merge_ref_alt
awk 'FNR==NR{a[$1]=$2FS$3FS$4FS$5FS$6FS$7FS$8FS$9FS$10FS$11FS$12FS$13FS$14FS$15FS$16FS$17FS$18; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, a[$1]}
' $path_ftemp_merge_alt_ref_clean $path_ftemp_dbsnp_extraction_ref_alt > $path_ftemp_merge_ref_alt

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 2:"
  head -10 $path_ftemp_merge_ref_alt
  echo "----------"
  echo "----------"
  echo "----------"
fi

# $1               $2 $3    $4  $5  $6  $7    $8    $9     $10    $11 $12 $13 $14 $15 $16  $17  $18 $19 $20 $21 $22  $23   $24
# identifier_merge ID CHROM POS ALT REF RS_ID ID_ar ALT_ar REF_ar SNP CHR BP  A1  A2  A1AF BETA SE  P   N   Z   INFO NCASE NCONT

# Check and filter the information from the merge.
echo "identifier_merge ID_ra ALT_ra REF_ra ID_ar ALT_ar REF_ar SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_merge_ref_alt_clean
#cat $path_ftemp_merge_ref_alt | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_ftemp_merge_ref_alt_clean
cat $path_ftemp_merge_ref_alt | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 24)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( $1 == "NA" )
    # Missing identifier indicates that second file did not match first file.
    # Skip the record.
    next
  else if (($12 != $3) || ($13 != $4) || ((toupper($14) != toupper($5)) && (toupper($14) != toupper($6))))
    # Chromosome, position, or alleles do not match.
    # Notice the comparison of A1 to both ALT and REF.
    next
  else
    # Keep record from successful merge.
    print $1, $2, $5, $6, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24
}' >> $path_ftemp_merge_ref_alt_clean

# $1               $2    $3     $4     $5    $6     $7     $8  $9  $10 $11 $12 $13  $14  $15 $16 $17 $18 $19  $20   $21
# identifier_merge ID_ra ALT_ra REF_ra ID_ar ALT_ar REF_ar SNP CHR BP  A1  A2  A1AF BETA SE  P   N   Z   INFO NCASE NCONT

# Determine whether to keep information from Merge 1 or Merge 2.
echo "identifier_merge ID ALT REF SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_merge_priority
cat $path_ftemp_merge_ref_alt_clean | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 21)
    # Skip any rows with incorrect count of column fields.
    next
  else if (($11 == $6) && ($12 == $7))
    # A1 = ALT
    # A2 = REF
    # Use Merge 1: ALT_REF
    print $1, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21
  else if (($11 == $4) && ($12 == $3))
    # A1 = REF
    # A2 = ALT
    # Use Merge 2: REF_ALT
    print $1, $2, $3, $4, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21
  else
    # Skip any row for which chromosome, position, or alleles do not match.
    next
}' >> $path_ftemp_merge_priority

# $1               $2 $3  $4  $5  $6  $7 $8 $9 $10  $11  $12 $13 $14 $15 $16  $17   $18
# identifier_merge ID ALT REF SNP CHR BP A1 A2 A1AF BETA SE  P   N   Z   INFO NCASE NCONT

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 1, Merge 2, and prioritization:"
  head -10 $path_ftemp_merge_priority
  echo "----------"
  echo "----------"
  echo "----------"
fi



##########
# 5. Adjust format of product GWAS summary statistics.

echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_product_format
cat $path_ftemp_merge_priority | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print $2, $6, $7, toupper($8), toupper($9), $10, $11, $12, $13, $14, $15, $16, $17, $18
}' >> $path_ftemp_product_format

# Compress file format.
gzip -cvf $path_ftemp_product_format > $path_file_gwas_product



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Filter and constrain values in GWAS summary statistics."
  echo "----------"
  echo "path to source GWAS file: " $path_file_gwas_source
  echo "path to product GWAS file: " $path_file_gwas_product
  echo "----------"
  echo "table before transformation:"
  zcat $path_file_gwas_source | head -5
  echo "- - Count of lines in source GWAS summary statistics:"
  zcat $path_file_gwas_source | wc -l
  echo "----------"
  echo "table after transformation:"
  zcat $path_file_gwas_product | head -5
  echo "- - Count of lines in product GWAS summary statistics:"
  zcat $path_file_gwas_product | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi



################################################################################
# Remove temporary, intermediate directories and files.

# Suppress this block for debugging.
if false; then
  rm -r $path_directory_temporary
fi



#
