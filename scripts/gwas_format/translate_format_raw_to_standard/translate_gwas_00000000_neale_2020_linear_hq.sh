#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Neale Lab, 2020.
# Host: https://pan.ukbb.broadinstitute.org/downloads

# Source Format
# Documentation: https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files#per-phenotype-files
# Human Genome Assembly: GRCh37 (UK Biobank)
# Effect Allele: "alt" ("Alternate allele ... Used as effect allele for GWAS.") (TCW; 16 February 2023)
# Delimiter: white space
# Columns: chr pos ref alt af_meta_hq beta_meta_hq se_meta_hq neglog10_pval_meta_hq neglog10_pval_heterogeneity_hq af_meta beta_meta se_meta neglog10_pval_meta neglog10_pval_heterogeneity ...
#          1   2   3   4   5          6            7          8                     9                              10      11        12      13                 14                          (TCW; 16 February 2023)
# Note: Use the general columns for Effect Allele Frequency (af_meta), Beta Coefficient (beta_meta),
# Standard Error (se_meta), and Probability (pval_meta) would be more inclusive as some
# phenotypes did not pass quality control.

# Format Translation
# The GWAS summary statistics do not include rs identifiers for SNP variants.
# Instead, define for SNP variants a special identifier format, "chr[chromosome]_[position]_[effect allele]".
# The LDSC Munge procedure might be able to interpret this identifier.
# columns: ("chr"$1"_"$2"_"$4), ...

# Product Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# Review: TCW; 21 February 2023

###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
fill_observations=${3} # logical binary indicator of whether to fill count of observations across all variants
observations=${4} # count of observations
fill_case_control=${5} # logical binary indicator of whether to fill counts of cases and controls across all variants
cases=${6} # count of cases
controls=${7} # count of controls
report=${8} # whether to print reports

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
# For conciseness, only support the conditions that are relevant.
if [ "$fill_observations" == "1" ] && [ "$fill_case_control" != "1" ]; then
  zcat $path_file_source | awk -v observations=$observations 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    if ((toupper($13) != "NA") && (($13 + 0) > 0))
      print ("chr"$1"_"$2"_"$4), $1, $2, toupper($4), toupper($3), $10, $11, $12, (10^-$13), (observations), "NA", (1.0), "NA", "NA"
    else
      print ("chr"$1"_"$2"_"$4), $1, $2, toupper($4), toupper($3), $10, $11, $12, "NA", (observations), "NA", (1.0), "NA", "NA"
  }' >> $path_file_temporary_format
fi


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



#
