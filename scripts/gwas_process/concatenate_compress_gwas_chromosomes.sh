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
pattern_gwas_report_file=${1} # string glob pattern by which to recognize PLINK2 GWAS report files
path_gwas_source_parent=${2} # full path to parent directory for GWAS across chromosomes
path_gwas_concatenation=${3} # full path to file for concatenation of GWAS across chromosomes
path_gwas_concatenation_compress=${4} # full path to file for concatenation of GWAS across chromosomes
report=${5} # whether to print reports

################################################################################
# Organize variables.

###########################################################################
# Execute procedure.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "path to source directory: " $path_gwas_source_parent
  echo "path to target file: " $path_gwas_concatenation_compress
fi

# Remove any previous versions of target files.
rm $path_gwas_concatenation
rm $path_gwas_concatenation_compress

#matches=("${path_source_chromosome}/${pattern_gwas_report_file}")
#path_source_file="${matches[0]}"

# Initialize table columns by file for chromosome 1.
# echo "#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P" > $path_gwas_concatenation
path_source_chromosome="$path_gwas_source_parent/chromosome_1"
matches=$(find "${path_source_chromosome}" -name "$pattern_gwas_report_file")
path_source_file=${matches[0]}
echo "source file for table column headers: " $path_source_file
cat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_gwas_concatenation

# Concatenate GWAS reports from selection chromosomes.
#for (( index=$chromosome_start; index<=$chromosome_end; index+=1 )); do
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "x")
for chromosome in "${chromosomes[@]}"; do
  path_source_chromosome="$path_gwas_source_parent/chromosome_${chromosome}"
  matches=$(find "${path_source_chromosome}" -name "$pattern_gwas_report_file")
  path_source_file=${matches[0]}
  echo "source file: " $path_source_file
  # Concatenate information from chromosome reports.
  #cat $path_source_file >> $path_gwas_concatentation
  #cat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12 }' >> $path_gwas_concatenation
  cat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 { print $0 }' >> $path_gwas_concatenation
done
# Compress file format.
gzip -cvf $path_gwas_concatenation > $path_gwas_concatenation_compress

###########################################################################
# Remove temporary files.
rm $path_gwas_concatenation
