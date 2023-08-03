#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution: 3 August 2023
# Review: TCW; 3 August 2023
################################################################################
# Note

# This script calculates the effective sample size (N-effective) for the counts
# of observations on each variant (SNP) in a logistic GWAS.
# References:
# 1. Zhou, Nature Neuroscience, 2020 (PubMed:32451486)
# 2. Walters, Nature Neuroscience, 2018 (PubMed:30482948)
# 3. Willer, Bioinformatics, 2010 (PubMed:20616382)
# 4. https://nealelab.github.io/UKBB_ldsc/viz_sampsize.html

# The calculation in this script assumes that the values in column "N" are the
# count of total observations (sum of cases and controls) for each variant (SNP)
# and that the values in column "NCASE" are the count of observations of cases
# for each variant (SNP).
# Hence: N = NCASE + NCONT
# The calculation also assumes that the count of total observations includes
# all (100%) of cases and that any discrepancies in count of total observations
# are indicative of differences in the count of cases.
# In fact the script ("pipe_gwas_clean.sh") after processing by GWAS2VCF does
# perform the calculation "NCONT = N - NCASE".

# Source Format
# Effect allele: "A1"
# Delimiter: white space
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P N  Z  INFO NCASE NCONT
#          1   2   3  4  5  6    7    8  9 10 11 12   13    14

# Product Format: Team Standard
# Effect allele: "A1"
# Delimiter: white space
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P N  Z  INFO NCASE NCONT
#          1   2   3  4  5  6    7    8  9 10 11 12   13    14

################################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics with GZip compression
report=${3} # whether to print reports

################################################################################
# Organize paths.

# Directories.
name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}" # hopefully unique

# Files.
path_file_temporary_calculation="${path_directory_product_temporary}/${name_base_file_product}_calculation.txt"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Initialize files.
rm $path_file_product



###########################################################################
# Execute procedure.



##########
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
# Calculate effective sample size (N-effective).
# References:
# 1. Zhou, Nature Neuroscience, 2020 (PubMed:32451486)
# 2. Walters, Nature Neuroscience, 2018 (PubMed:30482948)
# 3. Willer, Bioinformatics, 2010 (PubMed:20616382)
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_calculation
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  if (((toupper($10) != "NA") && (($10 + 0) > 0)) && ((toupper($13) != "NA") && (($13 + 0) > 0)))
    # Counts of observations and cases have non-missing values.
    # Calculate effective sample size of observations.
    # N = NCASE + NCONT
    # NCONT = N - NCASE
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, (4/((1/$13)+(1/($10-$13)))), $11, $12, $13, $14
  else
    # Counts of observations and cases have missing values.
    # Print the row entirely without change.
    print $0
}' >> $path_file_temporary_calculation



##########

# Compress file format.
gzip -cvf $path_file_temporary_calculation > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "calculate_effective_observations_logistic_gwas.sh"
  echo "----------"
  echo "Calculate effective sample size (N-effective)"
  echo "for summary statistics of logistic GWAS."
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "table after calculation:"
  head -10 $path_file_temporary_calculation
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
