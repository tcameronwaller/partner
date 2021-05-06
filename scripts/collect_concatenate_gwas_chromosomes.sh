#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize arguments.
pattern_source_file=${1} # glob pattern by which to recognize relevant files in source directory
path_source_directory=${2} # full path to parent directory for GWAS across chromosomes
chromosome_start=${3}
chromosome_end=${4}
path_gwas_concatenation=${5} # full path to file for concatenation of GWAS across chromosomes
path_gwas_concatenation_compress=${6} # full path to file for concatenation of GWAS across chromosomes
report=${7} # whether to print reports

################################################################################
# Organize variables.

###########################################################################
# Execute procedure.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "path to original directory: " $path_source_directory
  echo "path to new file: " $path_gwas_concatenation
  echo "start chromosome: " $chromosome_start
  echo "end chromosome: " $chromosome_end
fi

# Remove any previous versions of temporary files.
rm $path_gwas_concatenation

# Concatenate GWAS reports from selection chromosomes.
for (( index=$chromosome_start; index<=$chromosome_end; index+=1 )); do
  path_source_chromosome="$path_source_directory/chromosome_${index}"
  matches=("${path_source_chromosome}/${pattern_source_file}")
  path_source_file=${matches[0]}
  echo "source file: " $path_source_file
  # Concatenate information from chromosome reports.
  #cat $path_source_file >> $path_gwas_concatentation
done
# Compress file format.
#gzip -cvf $path_gwas_concatenation > $path_gwas_concatenation_compress

###########################################################################
# Remove temporary files.
#rm $path_gwas_concatenation
