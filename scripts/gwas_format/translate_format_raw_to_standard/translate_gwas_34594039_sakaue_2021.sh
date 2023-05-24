#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Notes

# This script translates the format of GWAS summary statistics from
# Sakaue et al, Nature Genetics, 2021 (PubMed:34594039).
# Host: https://www.ebi.ac.uk/gwas/publications/34594039

# Source Format
# Human Genome Assembly: GRCh37 (UK Biobank)
# Effect Allele: "effect_allele"
# Delimiter: white space
# Columns: chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error p_value variant_id
#          1          2                  3             4            5                       6    7              8       9          (TCW; 15 February 2023)

# Product Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT

# Review: TCW; 15 February 2023


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
if [ "$fill_observations" == "1" ] && [ "$fill_case_control" == "1" ]; then
  zcat $path_file_source | awk -v observations=$observations -v cases=$cases -v controls=$controls 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    print $9, $1, $2, toupper($3), toupper($4), $5, $6, $7, $8, (observations), "NA", (1.0), (cases), (controls)
  }' >> $path_file_temporary_format
fi

# Compress file format.
gzip -cvf $path_file_temporary_format > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Translate format of GWAS summary statistics."
  echo "----------"
  echo "path to source GWAS file: " $path_file_source
  echo "path to product GWAS file: " $path_file_product
  echo "----------"
  echo "table before format translation:"
  zcat $path_file_source | head -5
  echo "----------"
  echo "table after format translation:"
  zcat $path_file_product | head -5
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
