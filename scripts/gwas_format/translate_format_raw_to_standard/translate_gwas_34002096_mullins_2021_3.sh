#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 2 August 2023
# Date, last execution: 2 August 2023
# Review: TCW; 19 December 2023
################################################################################
# Note

################################################################################

# Maximal value in column "$14" (cases): 25,207
# Maximal value in column "$15" (controls): 515,930
# Maximal value in column "$13" (half effective sample size): 32,694.69

# Determine how many lines to skip.
# $ zcat ./pgc-bip2021-BDI.vcf.tsv.gz | awk 'BEGIN{FS=" "; OFS=" "} NR>72{print $0}' | head -10

# Determine maximal counts of observations (samples) for cases and controls.
# $ zcat pgc-bip2021-BDI.vcf.tsv.gz | awk 'BEGIN{FS=" "; OFS=" "; a=0} NR>73{if ((toupper($13) != "NA") && (($13+0)>(a+0))) a=$13} END{print a}'
# $ zcat pgc-bip2021-BDI.vcf.tsv.gz | awk 'BEGIN{FS=" "; OFS=" "; a=0} NR>73{if ((toupper($14) != "NA") && (($14+0)>(a+0))) a=$14} END{print a}'
# $ zcat pgc-bip2021-BDI.vcf.tsv.gz | awk 'BEGIN{FS=" "; OFS=" "; a=0} NR>73{if ((toupper($15) != "NA") && (($15+0)>(a+0))) a=$15} END{print a}'

# Count of cases: 25,060
# Count of controls: 449,978
# Count of total observations: 475,038
# Proportion of cases: 0.053
# Proportion of controls: 0.947

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

##########
# Translate format of GWAS summary statistics.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_format
# For conciseness, only support the conditions that are relevant.
if [ "$fill_observations" != "1" ] && [ "$fill_case_control" != "1" ]; then
  zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 73 {
    print $3, $1, $2, toupper($4), toupper($5), (($10*0.053)+($11*0.947)), $6, $7, $8, ($14 + $15), "NA", $12, $14, $15
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
