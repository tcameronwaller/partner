#!/bin/bash

###########################################################################
###########################################################################
###########################################################################

# "Organize logistic GWAS summary statistics in format for analysis by team."

# Source Format
# effect allele: "A1"

# Product Format (Team Standard)
# effect allele: "A1"
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# review: TCW; 10 August 2022

###########################################################################
###########################################################################
###########################################################################
# Notes

################################################################################
# Organize arguments.

path_file_gwas_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_gwas_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
report=${3} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_gwas_product .txt.gz)"
path_directory_product="$(dirname $path_file_gwas_product)"
path_directory_product_temporary="${path_directory_product}/temporary_format_${name_base_file_product}" # hopefully unique

path_file_temporary_constraint="${path_directory_product_temporary}/${name_base_file_product}_constraint.txt"
path_file_temporary_format="${path_directory_product_temporary}/${name_base_file_product}_format.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_gwas_product

###########################################################################
# Execute procedure.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "path to source GWAS file: " $path_file_gwas_source
  echo "path to product GWAS file: " $path_file_gwas_product
fi

# Here are some useful regular expressions to evaluate values in "awk".
# ( $12 ~ /^[0-9]+$/ ); ( $12 ~ /^[[:alpha:]]+$/ ); ( $12 ~ /^[[:punct:]]+$/ )

##########
# Constrain GWAS probability values within bounds for LDSC.
# Probability values in column "P" from PLINK2 have values of "NA" for missing
# or floating point values between zero and one.
# LDSC uses 64-bit floating point precision (double precision) to represent
# values from +/- 2.23E-308 to +/- 1.80E308 (https://github.com/bulik/ldsc/issues/144).
# Constrain probability values from 1.0E-305 to 1.0.
zcat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary_constraint
zcat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 15)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ( $15 != "NA" ) && ( ($15 + 0) < 1.0E-305 ) )
    # Constrain probability value.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, ( 1.0E-305 )
  else if ( ( $15 != "NA" ) && ( ($15 + 0) > 1.0 ) )
    # Constrain probability value.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, ( 1.0 )
  else
    print $0
  }' >> $path_file_temporary_constraint

##########
# Format of GWAS summary statistics for Team.

# Format of GWAS reports by PLINK2 for linear regression (".glm.linear").
# https://www.cog-genomics.org/plink/2.0/formats
# PLINK2 report format is similar for logistic regression (".glm.logistic").
# Logistic report has "OR" in place of "BETA", but positions are the same.

# PLINK2's format documentation describes column "REF" as "Reference allele" and
# column "ALT" as "All alternate alleles, comma-separated". Column "REF" and
# column "ALT" are always different (TCW, 7 July 2021).
# PLINK2's format documentation also describes column "A1" as "Counted allele in
# regression". Column "A1" always matches either column "REF" or column "ALT"
# (TCW, 7 July 2021).

# The designations of "reference" and "alternate" alleles are irrelevant to
# downstream analyses on GWAS summary statistics (such as LDSC).
# What matters is which allele PLINK2 counted as the "effect" allele in GWAS
# regression and which allele was the other allele.

# Column "A1" is the "effect" allele that PLINK2 counted in the regression,
# corresponding to the coefficient (beta).
# Determine the "non-effect" or "other" allele from either the "reference" or
# "alternate" allele, whichever differs from the "effect" allele.

# description: ............................ Team column ........... PLINK column .......... position
# variant identifier (RS ID): .............  "SNP" ................  "ID" ..................  3
# chromosome: .............................  "CHR" ................  "CHROM" ...............  1
# base pair coordinate: ...................  "BP" .................  "POS" .................  2
# alternate allele (effect allele): .......  "A1" .................  "A1" ..................  6
# reference allele (non-effect allele): ...  "A2" .................  "REF" or "ALT" ........  4 or 5
# alternate allele frequency: .............  "A1AF" ...............  "A1_FREQ" .............  7
# effect (coefficient): ...................  "BETA" ...............  "OR" .................. 12
# effect standard error: ..................  "SE" .................  "LOG(OR)_SE" .......... 13
# probability (p-value): ..................  "P" ..................  "P" ................... 15
# count total samples: ....................  "N" ..................  "OBS_CT" .............. 11
# probability (p-value): ..................  "Z" ..................  "Z_STAT" .............. 14
# imputation quality score: ...............  "INFO" ...............  "NA" .................. NA
# count case samples: .....................  "NCASE" ..............  "NA" .................. NA
# count control samples: ..................  "NCONT" ..............  "NA" .................. NA

# Organize information from linear GWAS.
# Select relevant columns and place them in the correct order.
# Convert original 'odds ratio' to 'log odds' or 'coefficient'.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
cat $path_file_temporary_constraint | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ($6 == $5 && $6 != $4)
    print $3, $1, $2, toupper($6), toupper($4), $7, (log($12)), $13, $15, $11, $14, (1), "NA", "NA"
  else if ($6 == $4 && $6 != $5)
    print $3, $1, $2, toupper($6), toupper($5), $7, (log($12)), $13, $15, $11, $14, (1), "NA", "NA"
  else
    next
  }' >> $path_file_temporary_format

# Compress file format.
gzip -cvf $path_file_temporary_format > $path_file_gwas_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "after bound constraints:"
  head -10 $path_file_temporary_constraint
  echo "after format:"
  head -10 $path_file_temporary_format
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
