#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 5 July 2024
# Date, last execution or modification: 9 July 2024
# Review: TCW; 9 July 2024
###############################################################################
# Note

# For accuracy in this file conversion from CRAM format to BAM format, the
# reference genome sequence must be the same as was used for the initial
# conversion from BAM format to CRAM format. Consult the header of the file in
# CRAM format to confirm which reference genome sequence was used.

################################################################################


################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source genomic or transcriptomic sequence information in CRAM format
path_file_product=${2} # full path to file for product genomic or transcriptomic sequence information in BAM format
path_file_reference_genome=${3} # full path to file for reference genome sequence
path_file_reference_genome_index=${4} # full path to file for reference genome sequence index
path_execution_samtools=${5} # full path to executable file for SamTools
report=${6} # whether to print reports

################################################################################
# Organize paths.

################################################################################
# Organize paths.

stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_product .bam)"
path_directory_product="$(dirname $path_file_product)"
#path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique
#path_file_temporary_1="${path_directory_product_temporary}/${name_base_file_product}_temporary_1.txt"

# Initialize directory.
mkdir -p $path_directory_product
#rm -r $path_directory_product_temporary
#mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###############################################################################
# Organize parameters.


###############################################################################
# Execute procedure.

$path_execution_samtools \
view \
-T $path_file_reference_genome \
-t $path_file_reference_genome_index \
--bam \
-o $path_file_product \
$path_file_source

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
#rm -r $path_directory_product_temporary

###############################################################################
# End.
