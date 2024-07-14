#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 July 2024
# Date, last execution or modification: 12 July 2024
# Review: TCW; 12 July 2024
###############################################################################
# Note

# http://www.htslib.org/doc/samtools-view.html
# http://www.htslib.org/doc/samtools-flags.html

###############################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source genomic or transcriptomic sequence information in BAM format
path_file_product=${2} # full path to file for product genomic or transcriptomic sequence information in BAM format
threads=${3} # count of concurrent or parallel process threads on node cores
report=${4} # whether to print reports to terminal
path_execution_samtools=${5} # full path to executable file for SamTools

###############################################################################
# Organize paths.

stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_product .bam)"
path_directory_product="$(dirname $path_file_product)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique
path_file_temporary_1="${path_directory_temporary}/${name_base_file_product}_temporary_1.bam"
path_file_product_index="${path_directory_product}/${name_base_file_product}.bam.bai"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_temporary

# Remove any previous version of the product file.
rm $path_file_temporary_1
rm $path_file_product
rm $path_file_product_index

###############################################################################
# Organize parameters.


###############################################################################
# Execute procedure.

if true; then
  # Filter.
  $path_execution_samtools \
  view \
  --bam \
  --with-header \
  --require-flags 0x1,0x2 \
  --excl-flags 0x4,0x8,0x200 \
  --threads $threads \
  --output $path_file_temporary_1 \
  $path_file_source
  # Sort coordinates of file in BAM format.
  $path_execution_samtools \
  sort \
  -n \
  -@ $threads \
  -o $path_file_product \
  $path_file_temporary_1
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
  echo "Filter reads, sort reads by name, and create index for file in BAM format."
  echo "----------"
  echo "path to source file: " $path_file_source
  echo "path to product sort file: " $path_file_product
  echo "path to product index file: " $path_file_product_index
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary

###############################################################################
# End.
