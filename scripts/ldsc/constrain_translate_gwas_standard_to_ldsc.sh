#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 23 December 2022
# Date, last execution: 29 November 2023
# Review: TCW; 29 November 2023
################################################################################
# Note

# Source Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# Product Format (LDSC)
# delimiter: white space
# columns: SNP A1 A2 N BETA P

################################################################################



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

path_file_temporary_constraint="${path_directory_product_temporary}/${name_base_file_product}_constraint.txt"
path_file_temporary_format="${path_directory_product_temporary}/${name_base_file_product}_format.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###########################################################################
# Execute procedure.

# Here are some useful regular expressions to evaluate values in "awk".
# ( $12 ~ /^[0-9]+$/ ); ( $12 ~ /^[[:alpha:]]+$/ ); ( $12 ~ /^[[:punct:]]+$/ )

##########
# Constrain GWAS probability values within bounds for LDSC.
# Probability values in column "P" from PLINK2 have values of "NA" for missing
# or floating point values between zero and one.
# Allow the LDSC Munge procedure to remove table rows for SNPs with missing
# values ("NA") in columns for beta coefficient, standard error, or probability.
# LDSC uses 64-bit floating point precision (double precision) to represent
# values from +/- 2.23E-308 to +/- 1.80E308 (https://github.com/bulik/ldsc/issues/144).
# Constrain probability values from 1.0E-305 to 1.0.

zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary_constraint
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 14)
    # Skip any of table rows with incorrect count of column fields, indicating empty cells.
    next
  else if ( ($9 == "") || (toupper($9) == "NA") || (toupper($9) == "NAN") || ( ($9 + 0) < 0 ) )
    # Skip any rows with missing or nonsense probability.
    next
  else if ( ( ($9 + 0) > 0 ) && ( ($9 + 0) < 1.0E-307 ) )
    # Constrain probability.
    print $1, $2, $3, $4, $5, $6, $7, $8, (1.0E-307), $10, $11, $12, $13, $14
  else if ( ( ($9 + 0) > 1.0 ) )
    # Constrain probability.
    print $1, $2, $3, $4, $5, $6, $7, $8, (1.0), $10, $11, $12, $13, $14
  else
    # Row passes all checks, filters, and constraints.
    # Print the row entirely.
    print $0
  }' >> $path_file_temporary_constraint


##########
# Translate format of GWAS summary statistics.
# AWK interprets a single space delimiter (FS=" ") as any white space.
# LDSC seems to accommodate unique identifiers for Single Nucleotide
# Polymorphisms that include "rsIDs" or unique identifiers composed of the
# chromosome, base pair position, and allele (such as
# "chromosome_position_allele").
# effect allele: "A1"
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
# https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format

echo "SNP A1 A2 N BETA P" > $path_file_temporary_format
cat $path_file_temporary_constraint | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  print $1, toupper($4), toupper($5), $10, $7, $9
}' >> $path_file_temporary_format

# Compress file format.
gzip -cvf $path_file_temporary_format > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate format of GWAS summary statistics for LDSC."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "table after constraints:"
  head -10 $path_file_temporary_constraint
  echo "table after format:"
  head -10 $path_file_temporary_format
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
