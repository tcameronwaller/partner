#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 July 2024
# Date, last execution or modification: 11 July 2024
# Review: TCW; 11 July 2024
###############################################################################
# Note


###############################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source genomic or transcriptomic sequence information in CRAM format
path_file_product=${2} # full path to file for product genomic or transcriptomic sequence information in BAM format
threads=${3} # count of concurrent or parallel process threads on node cores
report=${4} # whether to print reports to terminal
path_execution_samtools=${5} # full path to executable file for SamTools

###############################################################################
# Organize paths.

#stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_product .bam)"
path_directory_product="$(dirname $path_file_product)"
path_file_product_index="${path_directory_product}/${name_base_file_product}.bam.bai"

# Initialize directory.
mkdir -p $path_directory_product

# Remove any previous version of the product file.
rm $path_file_product
rm $path_file_product_index

###############################################################################
# Organize parameters.


###############################################################################
# Execute procedure.

if true; then
  # Sort coordinates of file in BAM format.
  $path_execution_samtools \
  sort \
  -@ $threads \
  -o $path_file_product \
  $path_file_source
  # Create index for file in BAM format.
  $path_execution_samtools \
  index \
  --bai \
  --threads $threads \
  -o $path_file_product_index \
  $path_file_product
fi


##########
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Sort coordinates and create index for file in BAM format."
  echo "----------"
  echo "path to source file: " $path_file_source
  echo "path to product sort file: " $path_file_product
  echo "path to product index file: " $path_file_product_index
  echo "----------"
fi

###############################################################################
# End.
