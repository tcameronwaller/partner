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
  echo "path to new file: " $path_gwas_concatenation_compress
  echo "start chromosome: " $chromosome_start
  echo "end chromosome: " $chromosome_end
fi

# Remove any previous versions of temporary files.
rm $path_gwas_concatenation
rm $path_gwas_concatenation_compress

# TODO: make this more versatile by reading in the first row of the first chromosome file...

# Concatenate GWAS reports from selection chromosomes.
echo "#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P" > $path_gwas_concatenation
for (( index=$chromosome_start; index<=$chromosome_end; index+=1 )); do
  path_source_chromosome="$path_source_directory/chromosome_${index}"
  matches=("${path_source_chromosome}/${pattern_source_file}")
  path_source_file="${matches[0]}"
  echo "source file: " $path_source_file
  # Concatenate information from chromosome reports.
  #cat $path_source_file >> $path_gwas_concatentation
  cat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >> $path_gwas_concatenation
done
# Compress file format.
gzip -cvf $path_gwas_concatenation > $path_gwas_concatenation_compress

###########################################################################
# Remove temporary files.
rm $path_gwas_concatenation
