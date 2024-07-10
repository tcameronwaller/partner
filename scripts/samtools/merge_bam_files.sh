#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 10 July 2024
# Date, last execution or modification: 10 July 2024
# Review: TCW; 10 July 2024
###############################################################################
# Note

# The file in "path_file_source" is in text format and gives a list of paths to
# the multiple source files in BAM format to merge together.

###############################################################################
# Organize arguments.

path_file_source=${1} # full path to file in text format with list of paths to files in BAM format for merge
path_file_product=${2} # full path to file for product genomic or transcriptomic sequence information in BAM format
threads=${3} # count of concurrent or parallel process threads on node cores
report=${4} # whether to print reports to terminal
path_execution_samtools=${5} # full path to executable file for SamTools

###############################################################################
# Organize paths.

stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_product .bam)"
path_directory_product="$(dirname $path_file_product)"
#path_directory_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique
#path_file_temporary_1="${path_directory_temporary}/${name_base_file_product}_temporary_1.bam"

# Initialize directory.
mkdir -p $path_directory_product
#rm -r $path_directory_temporary
#mkdir -p $path_directory_temporary

# Remove any previous version of the product file.
#rm $path_file_temporary_1
rm $path_file_product

###############################################################################
# Organize parameters.


###############################################################################
# Execute procedure.

if true; then
  # Convert file from CRAM format to BAM format.
  $path_execution_samtools \
  merge \
  --threads $threads \
  -b $path_file_source \
  -o $path_file_product
fi


##########
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Merge together multiple files in BAM format."
  echo "----------"
  echo "path to source file: " $path_file_source
  echo "path to product file: " $path_file_product
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
#rm -r $path_directory_temporary

###############################################################################
# End.
