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
pattern_gwas_report_file=${1} # string glob pattern by which to recognize PLINK2 GWAS summary statistics report files
name_report_log_file=${2} # string name of file for PLINK2 report log
name_report_allele_frequency_file=${3} # string name of file for PLINK2 report on allele frequencies
path_gwas_source_parent=${4} # full path to parent directory for GWAS across chromosomes
path_gwas_target_parent=${5} # full path to parent directory for target GWAS summary statistics
path_gwas_concatenation=${6} # full path to file for concatenation of GWAS across chromosomes
path_gwas_concatenation_compress=${7} # full path to file for concatenation of GWAS across chromosomes
chromosome_x=${8} # whether to collect GWAS summary statistics report for Chromosome X
report=${9} # whether to print reports

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

# Define name of target file for PLINK2 Allele Frequency report.
path_target_file_frequency="${path_gwas_target_parent}/${name_report_allele_frequency_file}"

# Remove any previous versions of target files.
rm $path_target_file_frequency
rm $path_gwas_concatenation
rm $path_gwas_concatenation_compress

# Initialize table columns of PLINK2 Allele Frequencies by file for chromosome 1.
path_source_chromosome="$path_gwas_source_parent/chromosome_1"
path_source_file_frequency="${path_source_chromosome}/${name_report_allele_frequency_file}"
echo "source Allele Frequency file for table column headers: " $path_source_file_frequency
cat $path_source_file_frequency | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > ${path_target_file_frequency}

# Initialize table columns of GWAS summary statistics by file for chromosome 1.
# echo "#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P" > $path_gwas_concatenation
path_source_chromosome="$path_gwas_source_parent/chromosome_1"
#matches=("${path_source_chromosome}/${pattern_gwas_report_file}")
matches=$(find "${path_source_chromosome}" -name "$pattern_gwas_report_file")
path_source_file_gwas=${matches[0]}
echo "source GWAS Summary Statistics file for table column headers: " $path_source_file_gwas
cat $path_source_file_gwas | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_gwas_concatenation

# Concatenate GWAS reports from selection chromosomes.
#for (( index=$chromosome_start; index<=$chromosome_end; index+=1 )); do
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "x" "xy")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
for chromosome in "${chromosomes[@]}"; do
  # Define path to chromosome.
  path_source_chromosome="$path_gwas_source_parent/chromosome_${chromosome}"

  # Define source file for Plink2 report log.
  path_source_file_log="${path_source_chromosome}/${name_report_log_file}"
  # Copy file for Plink2 report log.
  path_target_file_log="${path_gwas_target_parent}/chromosome_${chromosome}_${name_report_log_file}"
  cp "${path_source_file_log}" "${path_target_file_log}"

  # Define source file for Allele Frequencies.
  path_source_file_frequency="${path_source_chromosome}/${name_report_allele_frequency_file}"
  # Concatenate information from chromosome reports on Allele Frequencies.
  cat $path_source_file_frequency | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 { print $0 }' >> ${path_target_file_frequency}

  # Define source file for GWAS Summary Statistics.
  matches=$(find "${path_source_chromosome}" -name "$pattern_gwas_report_file")
  path_source_file_gwas=${matches[0]}
  echo "GWAS source file: " $path_source_file_gwas
  # Concatenate information from chromosome reports on GWAS Summary Statistics.
  #cat $path_source_file >> $path_gwas_concatentation
  #cat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12 }' >> $path_gwas_concatenation
  cat $path_source_file_gwas | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 { print $0 }' >> $path_gwas_concatenation
done
# Compress file format.
gzip -vf $path_target_file_frequency
gzip -cvf $path_gwas_concatenation > $path_gwas_concatenation_compress

###########################################################################
# Remove temporary files.
rm $path_gwas_concatenation
