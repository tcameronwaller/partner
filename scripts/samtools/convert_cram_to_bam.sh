#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 5 July 2024
# Date, last execution or modification: 10 July 2024
# Review: TCW; 10 July 2024
###############################################################################
# Note

# For accuracy in this file conversion from CRAM format to BAM format, the
# reference genome sequence must be the same as was used for the initial
# conversion from BAM format to CRAM format. Consult the header of the file in
# CRAM format to confirm which reference genome sequence was used.


###############################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source genomic or transcriptomic sequence information in CRAM format
path_file_product=${2} # full path to file for product genomic or transcriptomic sequence information in BAM format
path_file_reference_genome=${3} # full path to file for sequence of reference genome
path_file_reference_genome_index=${4} # full path to index file for sequence of reference genome
threads=${5} # count of concurrent or parallel process threads on node cores
report=${6} # whether to print reports to terminal
path_execution_samtools=${7} # full path to executable file for SamTools

###############################################################################
# Organize paths.

stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_product .bam)"
path_directory_product="$(dirname $path_file_product)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique
path_file_temporary_1="${path_directory_temporary}/${name_base_file_product}_temporary_1.bam"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_temporary

# Remove any previous version of the product file.
rm $path_file_temporary_1
rm $path_file_product

###############################################################################
# Organize parameters.


###############################################################################
# Execute procedure.

if false; then
  # Convert file from CRAM format to BAM format.
  $path_execution_samtools \
  view \
  --bam \
  --with-header \
  --threads $threads \
  -T $path_file_reference_genome \
  -t $path_file_reference_genome_index \
  -o $path_file_product \
  $path_file_source
fi

if true; then
  # Convert file from CRAM format to BAM format.
  $path_execution_samtools \
  view \
  --bam \
  --with-header \
  --threads $threads \
  -T $path_file_reference_genome \
  -t $path_file_reference_genome_index \
  --output $path_file_temporary_1 \
  $path_file_source
  # Sort coordinates of file in BAM format.
  $path_execution_samtools \
  sort \
  -@ $threads \
  -o $path_file_product \
  $path_file_temporary_1
fi


##########
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Convert genomic or transcriptomic sequence data from CRAM to BAM file format."
  echo "----------"
  echo "path to source file: " $path_file_source
  echo "path to product file: " $path_file_product
  echo "----------"
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary

###############################################################################
# End.
