#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 1 March 2023
# Date, last execution: 18 December 2023
# Review: TCW; 18 December 2023
################################################################################
# Note

# This script translates the format of GWAS summary statistics from
# Walters et al, Nature Neuroscience, 2018
# (PubMed:30482948).
# Host: https://pgc.unc.edu/
# Host: https://figshare.com/articles/dataset/sud2018-alc/14672187


# pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt.gz
# Determine maximal counts of observations (samples) for cases and controls.
# $ zcat pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt.gz | awk 'BEGIN{FS=" "; OFS=" "; a=0} NR>1{if ((toupper($10) != "NA") && (($10+0)>(a+0))) a=$10} END{print a}'
# Maximal value in column "$10" (effective sample size): 28,664.4


################################################################################



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

################################################################################
# Execute procedure.

# Note: TCW; 17 February 2023
# The logarithm of a negative number or zero is undefined.

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
# For conciseness, only support the conditions that are relevant.
if [ "$fill_observations" == "1" ] && [ "$fill_case_control" == "1" ]; then
  zcat $path_file_source | awk -v observations=$observations -v cases=$cases -v controls=$controls 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    if ((toupper($7) != "NA") && (($7 + 0) > 0))
      print $2, $1, $3, toupper($4), toupper($5), "NA", log($7), $8, $9, (observations), "NA", $6, (cases), (controls)
    else
      print $2, $1, $3, toupper($4), toupper($5), "NA", "NA", $8, $9, (observations), "NA", $6, (cases), (controls)
  }' >> $path_file_temporary_format
elif [ "$fill_observations" != "1" ] && [ "$fill_case_control" == "1" ]; then
  zcat $path_file_source | awk -v cases=$cases -v controls=$controls 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    if ((toupper($7) != "NA") && (($7 + 0) > 0))
      print $2, $1, $3, toupper($4), toupper($5), "NA", log($7), $8, $9, int($10), "NA", $6, (cases), (controls)
    else
      print $2, $1, $3, toupper($4), toupper($5), "NA", "NA", $8, $9, int($10), "NA", $6, (cases), (controls)
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
