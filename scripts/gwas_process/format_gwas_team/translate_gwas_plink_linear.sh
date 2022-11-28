#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# Source Format
# delimiter: white space
# columns: CHROM POS ID REF ALT A1 A1_FREQ TEST OBS_CT BETA SE T_STAT P

# Format Translation
# columns: $3, $1, $2, $6, ($4 or $5), $7, $10, $11, $13, $9, "NA", (1), "NA", "NA"

# Product Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# review: TCW; 23 November 2022

###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
report=${3} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_format_${name_base_file_product}" # hopefully unique
path_file_temporary_format="${path_directory_product_temporary}/${name_base_file_product}_format.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###########################################################################
# Execute procedure.

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.

echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ($6 == $5 && $6 != $4)
    print $3, $1, $2, toupper($6), toupper($4), $7, $10, $11, $13, $9, "NA", (1), "NA", "NA"
  else if ($6 == $4 && $6 != $5)
    print $3, $1, $2, toupper($6), toupper($5), $7, $10, $11, $13, $9, "NA", (1), "NA", "NA"
  else
    next
  }' >> $path_file_temporary_format

# Compress file format.
gzip -cvf $path_file_temporary_format > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate format of GWAS summary statistics."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "table after format:"
  head -10 $path_file_temporary_format
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



##########
# Notes on format of GWAS summary statistics from PLINK2.

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
# effect (coefficient): ...................  "BETA" ...............  "BETA" ................ 10
# effect standard error: ..................  "SE" .................  "SE" .................. 11
# probability (p-value): ..................  "P" ..................  "P" ................... 13
# count total samples: ....................  "N" ..................  "OBS_CT" ..............  9
# probability (p-value): ..................  "Z" ..................  "NA" .................. (BETA / SE) # division by zero
# imputation quality score: ...............  "INFO" ...............  "NA" .................. NA
# count case samples: .....................  "NCASE" ..............  "NA" .................. NA
# count control samples: ..................  "NCONT" ..............  "NA" .................. NA



#
