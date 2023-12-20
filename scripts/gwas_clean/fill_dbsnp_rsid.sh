#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 December 2023
# Date, last execution: 19 December 2023
# Review: 19 December 2023
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

# In situations for which dbSNP references multiple rsIDs for a single SNP, this
# procedure keeps only the first rsID from the list.

# The script below is a companion to perform the first few operations of this
# script to save time in subsequent iterations.
# "/.../partner/scripts/bcftools/extract_dbsnp_biallelic_sites_allele_identifiers.sh"

# Note: TCW; 19 December 2023
# The process of the full "fill_dbsnp_rsid.sh" script, without any division into
# smaller sections (subscripts), completed for the study
# "37872160_williams_2023", using 278.20 Gigabytes of memory with management by
# the SLURM scheduler.
# Subsequently, the process of the "fill_dbsnp_rsid.sh" script, after division
# into five smaller sections (subscripts) in an attempt to save memory,
# completed for the study "37872160_williams_2023", using 274.89 Gigabytes of
# memory with management by the SLURM scheduler.
# It appears that division of the process into smaller subscripts did not offer
# any advantage of clearing system memory between subscripts.

# TODO: TCW; 13 December 2023
# Write out two different "product files"...
# 1. Filtered to include only SNPs that had matching rsIDs from dbSNP
#   - This table is useful for a count of the "successful" SNPs if nothing else.
# 2. All SNPs from the original GWAS sum stats regardless.
#   - Subsequent procedures might be able to use more SNPs.
# IDEA: --> introduce an argument "strict" to determine whether to write out all SNPs or only matching


##########
# Standard format of GWAS summary statistics.
# This is the obligatory format of the source and product GWAS summary
# statistics.
# File name suffix: ".txt.gz"
# File compression: Gzip
# Delimiter: white space (" ")
# Columns and sequence:
# $1  $2  $3 $4 $5 $6   $7   $8 $9 $10 $11 $12  $13   $14
# SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT

##########
# Limitations.
# 1. The current implementation of this script makes internal reference to
# anonymous hard-coded file paths for the sake of convenience when executing in
# a specific computational environment. The script could of course be adapted to
# accept these file path variables as parameters (arguments).
# 2. For simplicity, this procedure filters to biallelic sites before extracting
# information about SNPs and their rsIDs from the dbSNP reference.

##########
# Testing.
# TCW; 13 December 2023
# Version of dbSNP: build 155
# Accession of dbSNP: TCW; 2023-02-06
# Translation of chromosome designations in dbSNP: TCW; 2023-02-15
# I tested the procedure using the GWAS summary statistics for levels of
# testosterone in females from Ruth et al, 2020 (PubMed:32042192). This set of
# GWAS summary statistics already included rsIDs for most SNPs, and 99.765% of
# the rsIDs matched from dbSNP were identical to the original rsIDs. Hence this
# procedure has high accuracy in its matching of rsIDs from dbSNP. The main
# limitation of this procedure is the accommodation of biallelic sites only from
# dbSNP, and the test GWAS summary statistics had lost about half (47%) of their
# SNPs after this procedure. So in summary, this procedure matches rsIDs from
# dbSNP with high accuracy but also loses a considerable proportion of SNPs,
# many of which might be for multiallelic sites.


################################################################################


################################################################################
# Organize arguments.

path_file_gwas_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_gwas_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
strict=${3} # whether to return GWAS summary statistics filtered to SNPs with successful match to dbSNP rsID
report=${4} # whether to print reports


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
path_ftemp_merge_priority_clean="${path_directory_temporary}/merge_priority_clean.txt"
path_ftemp_merge_priority_clean_strict="${path_directory_temporary}/merge_priority_clean_strict.txt"

path_ftemp_merge_priority_check="${path_directory_temporary}/merge_priority_check.txt"
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

##########################
#######################
###################
################
#############
# Begin part 1

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

# End part 1
############
##############
#################
#######################
###############################


##########################
#######################
###################
################
#############
# Begin part 2



##########
# 3. For GWAS summary statistics, assemble information with unique identifiers
# specific to site (chromosome, position) and both effect and other alleles.
# Identifier format: <CHROM>_<POS>_<ALT>_<REF>
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
# These merges treat Table 2 as the priority table.
# The merge of Table 1 with Table 2 will include all records from Table 2 but
# only those records from Table 1 that have matching identifiers.
# These merges print each row of Table 2 regardless of whether there is a record
# with matching identifier in Table 1.

# Process explanation:
# 1. GNU awk reads lines from the first file (FNR==NR) into an array in memory
# (array "a") that uses values from the first column of the first file (column
# "identifier_merge") as a hashable index (a[$1]=...).
# 2. GNU awk then proceeds to the second file and prints line by line while
# including any of the hashable lines from the first file that have a matching
# index, this time from the first column of the second file (column
# "identifier_merge").

# Reference:
# https://stackoverflow.com/questions/32481877/what-are-nr-and-fnr-and-what-does-nr-fnr-imply


##########
# 4.1. Merge 1.
# Table 1: dbSNP extraction with identifier format <CHROM>_<POS>_<ALT>_<REF>.
# Table 1: 7 total columns with merge identifier in column 1.
# Table 2: GWAS summary statistics
# Table 2: 15 total columns with merge identifier in column 1.
# Merged table: 21 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF RS_ID" > $path_ftemp_merge_alt_ref
awk 'FNR==NR{a[$1]=$2FS$3FS$4FS$5FS$6FS$7; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, a[$1]
} END {
  delete a
}' $path_ftemp_dbsnp_extraction_alt_ref $path_ftemp_gwas_identifier_a1_a2 > $path_ftemp_merge_alt_ref

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "This is your chance to check merge accuracy!"
  echo "Table immediately after Merge 1:"
  head -10 $path_ftemp_merge_alt_ref
  cat $path_ftemp_merge_alt_ref | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi

# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16 $17   $18 $19 $20  $21
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID  CHROM POS ALT REF RS_ID

# Check and filter the information from the merge.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar" > $path_ftemp_merge_alt_ref_clean
#cat $path_ftemp_merge_alt_ref | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_ftemp_merge_alt_ref_clean
cat $path_ftemp_merge_alt_ref | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 21)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ($1 == "NA") )
    # Missing identifier indicates that records did not match or other problem.
    # Skip the record.
    next
  else if ( ($3 == "NA") && ($4 == "NA") && ($5 == "NA") && ($6 == "NA") )
    # Missingness of match criteria from Table 2 should not occur, but this
    # would indicate that the merge procedure had kept records from Table 1
    # that did not match records in Table 2.
    # Skip the record.
    next
  else
    # Keep record from successful merge.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20
}' >> $path_ftemp_merge_alt_ref_clean

##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 1 and subsequent filters:"
  head -10 $path_ftemp_merge_alt_ref_clean
  cat $path_ftemp_merge_alt_ref_clean | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi


# End part 2
############
##############
#################
#######################
###############################



##########################
#######################
###################
################
#############
# Begin part 3



##########
# 4.2. Merge 2.
# Table 1: dbSNP extraction with identifier format <CHROM>_<POS>_<REF>_<ALT>.
# Table 1: 7 total columns with merge identifier in column 1.
# Table 2: Table from Merge 1.
# Table 2: 20 total columns with merge identifier in column 1.
# Merged table: 26 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID CHROM POS ALT REF RS_ID" > $path_ftemp_merge_ref_alt
awk 'FNR==NR{a[$1]=$2FS$3FS$4FS$5FS$6FS$7; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, a[$1]
} END {
  delete a
}' $path_ftemp_dbsnp_extraction_ref_alt $path_ftemp_merge_alt_ref_clean > $path_ftemp_merge_ref_alt



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "This is your chance to check merge accuracy!"
  echo "Table immediately after Merge 2:"
  head -10 $path_ftemp_merge_ref_alt
  cat $path_ftemp_merge_ref_alt | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi

# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16   $17      $18    $19    $20    $21 $22   $23 $24 $25 $26
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID  CHROM POS ALT REF RS_ID

# Check and filter the information from the merge.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID_ra CHROM_ra POS_ra ALT_ra REF_ra" > $path_ftemp_merge_ref_alt_clean
#cat $path_ftemp_merge_ref_alt | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_ftemp_merge_ref_alt_clean
cat $path_ftemp_merge_ref_alt | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 26)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ($1 == "NA") )
    # Missing identifier indicates that records did not match or other problem.
    # Skip the record.
    next
  else if ( ($3 == "NA") && ($4 == "NA") && ($5 == "NA") && ($6 == "NA") )
    # Missingness of match criteria from Table 2 should not occur, but this
    # would indicate that the merge procedure had kept records from Table 1
    # that did not match records in Table 2.
    # Skip the record.
    next
  else
    # Keep record from successful merge.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25
}' >> $path_ftemp_merge_ref_alt_clean



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 2 and subsequent filters:"
  head -10 $path_ftemp_merge_ref_alt_clean
  cat $path_ftemp_merge_ref_alt_clean | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi


# End part 3
############
##############
#################
#######################
###############################


##########################
#######################
###################
################
#############
# Begin part 4




# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16   $17      $18    $19    $20    $21   $22      $23    $24    $25
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID_ar CHROM_ar POS_ar ALT_ar REF_ar ID_ra CHROM_ra POS_ra ALT_ra REF_ra



# Note: TCW; 12 December 2023
# Where there are multiple rsIDs for a SNP, dbSNP records use semicolon ";" delimited lists.
# Split by semicolon and take the first instance in the array to accommodate.

# Determine whether to keep information from Merge 1 or Merge 2.
# First print conditional block: A1 is ALT, A2 is REF, use Merge 1
# Second print conditional block: A1 is REF, A2 is ALT, use Merge 2
# https://www.gnu.org/software/gawk/manual/html_node/String-Functions.html
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority
cat $path_ftemp_merge_ref_alt_clean | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 25)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ($3 == $17) && ($4 == $18) && ($5 == $19) && ($6 == $20) )
    # A1 is ALT
    # A2 is REF
    # Use Merge 1: ALT_REF
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20
  else if ( ($3 == $22) && ($4 == $23) && ($5 == $25) && ($6 == $24) )
    # A1 is REF
    # A2 is ALT
    # Use Merge 2: REF_ALT
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $21, $22, $23, $24, $25
  else
    # The record did not match or merge with information from dbSNP.
    # Fill with missing information.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "NA", "NA", "NA", "NA", "NA"
}' >> $path_ftemp_merge_priority
# GNU Awk can only handle a few operations at a time.
# Split delimited lists of identifiers and keep only the first.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority_clean
cat $path_ftemp_merge_priority | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  (a = $16); split(a, b, ";"); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, b[1], $17, $18, $19, $20
}' >> $path_ftemp_merge_priority_clean

# Strict.
echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority_clean_strict
cat $path_ftemp_merge_priority_clean | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( ($16 == "NA") && ($17 == "NA") && ($18 == "NA") && ($19 == "NA") && ($20 == "NA") )
    # Skip any row records with missing information from merge.
    next
  else
    # Keep the row record.
    print $0
}' >> $path_ftemp_merge_priority_clean_strict


# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16 $17   $18 $19 $20
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID  CHROM POS ALT REF



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Table after Merge 1, Merge 2, and prioritization:"
  head -10 $path_ftemp_merge_priority_clean
  cat $path_ftemp_merge_priority_clean | wc -l
  echo "----------"
  echo "----------"
  echo "----------"
fi



##########
# 5. Adjust format of product GWAS summary statistics.

if [ "$strict" == "true" ]; then

  # The strict version has already been filtered to non-missing matches with
  # dbSNP.

  echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_product_format
  cat $path_ftemp_merge_priority_clean_strict | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    print $16, $3, $4, toupper($5), toupper($6), $7, $8, $9, $10, $11, $12, $13, $14, $15
  }' >> $path_ftemp_product_format

elif [ "$strict" == "false" ]; then

  # The non-strict version has not yet been filtered to non-missing matches with
  # dbSNP.
  # Only replace the original SNP identifier if the match from dbSNP is not
  # missing.

    echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_product_format
    cat $path_ftemp_merge_priority_clean | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
      if ( ($16 != "NA") && ($17 != "NA") && ($18 != "NA") && ($19 != "NA") && ($20 != "NA") )
        print $16, $3, $4, toupper($5), toupper($6), $7, $8, $9, $10, $11, $12, $13, $14, $15
      else
        print $2, $3, $4, toupper($5), toupper($6), $7, $8, $9, $10, $11, $12, $13, $14, $15
    }' >> $path_ftemp_product_format

fi

# Compress file format.
gzip -cvf $path_ftemp_product_format > $path_file_gwas_product

# End part 4
############
##############
#################
#######################
###############################




##########################
#######################
###################
################
#############
# Begin part 5



# $1               $2  $3  $4 $5 $6 $7   $8   $9 $10 $11 $12 $13  $14   $15   $16 $17   $18 $19 $20
# identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT ID  CHROM POS ALT REF



if true; then
  # Test the procedure by determining the proportion of SNP rsIDs from dbSNP that match those in original GWAS summary statistics.
  echo "identifier_merge SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT ID CHROM POS ALT REF" > $path_ftemp_merge_priority_check
  cat $path_ftemp_merge_priority_clean_strict | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
    if ( $2 != $16)
      # Keep rows for which the source SNP rsID does not match the dbSNP rsID.
      print $0
    else
      # Skip rows for which source SNP rsID matches the dbSNP rsID.
      next
  }' >> $path_ftemp_merge_priority_check
  # Report.
  if [[ "$report" == "true" ]]; then
    # Calculations.
    # https://www.gnu.org/software/bc/manual/html_mono/bc.html
    count_check=$(cat $path_ftemp_merge_priority_check | wc -l)
    count_total=$(cat $path_ftemp_merge_priority_clean_strict | wc -l)
    proportion_check=$(echo "scale=5; ($count_check / $count_total)" | bc -l) # order of operations and scale rounding matters
    percentage_check=$(echo "scale=5; ($proportion_check * 100)" | bc -l) # order of operations and scale rounding matters
    echo "----------"
    echo "----------"
    echo "----------"
    echo "For special situations where the source GWAS summary statistics"
    echo "already included rsIDs for SNPs, it is possible to test the"
    echo "proportion of the rsIDs extracted from dbSNP that match the"
    echo "originals."
    echo "----------"
    echo "This check uses the strict version that only considers the SNPs that"
    echo "matched and merged with dbSNP rsIDs successfully."
    echo "----------"
    echo "Table that only includes SNPs for which the dbSNP rsID does not match"
    echo "the source rsID:"
    echo "----------"
    head -5 $path_ftemp_merge_priority_check
    echo "----------"
    echo "Lines that do not match between original and dbSNP: " $count_check
    echo "Lines total (strict): " $count_total
    echo "Proportion that do not match: " $proportion_check
    echo "Percentage that do not match: " $percentage_check "%"
    echo "----------"
    echo "----------"
    echo "----------"
  fi
fi





################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  # Calculations.
  # https://www.gnu.org/software/bc/manual/html_mono/bc.html
  count_source=$(zcat $path_file_gwas_source | wc -l)
  count_strict=$(cat $path_ftemp_merge_priority_clean_strict | wc -l)
  count_product=$(zcat $path_file_gwas_product | wc -l)
  #percentage_poduct=$(echo "100 * ($count_product / $count_source)" | bc -l) # order of operations and scale rounding matters
  proportion_strict=$(echo "scale=5; ($count_strict / $count_source)" | bc -l) # order of operations and scale rounding matters
  percentage_strict=$(echo "scale=5; ($proportion_strict * 100)" | bc -l) # order of operations and scale rounding matters
  count_difference=$(($count_source - $count_strict))
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Fill dbSNP rsIDs for SNPs in GWAS summary statistics."
  echo "----------"
  echo "Count of original lines in source: " $count_source
  echo "Count of lines that matched and merged with dbSNP (strict): " $count_strict
  echo "Count of lines difference between source and strict: " $count_difference
  echo "Proportion of lines that matched and merged: " $proportion_strict
  echo "Percentage of lines that matched and merged: " $percentage_strict "%"
  echo "----------"
  echo "Use strict filter to SNPs that matched and merged with dbSNP: " $strict
  echo "If active this filter determines count of lines in product."
  echo "----------"
  echo "path to source GWAS file: " $path_file_gwas_source
  echo "path to product GWAS file: " $path_file_gwas_product
  echo "----------"
  echo "source table before transformation:"
  zcat $path_file_gwas_source | head -5
  echo "- - Count of lines in source table: " $count_source
  echo "----------"
  echo "product table after transformation:"
  zcat $path_file_gwas_product | head -5
  echo "- - Count of lines in product table: " $count_product
  echo "----------"
  echo "----------"
  echo "----------"
fi



################################################################################
# Remove temporary, intermediate directories and files.

# Suppress this block for debugging.
if true; then
  rm -r $path_directory_temporary
fi



#
