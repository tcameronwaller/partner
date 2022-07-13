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
pattern_file_gwas_source=${1} # string glob pattern by which to recognize GWAS summary statistics report files from PLINK2
pattern_file_frequency_source=${2} # string glob pattern by which to recognize allele frequency report files from PLINK2
pattern_file_log_source=${3} # string glob pattern by which to recognize execution log files from PLINK2
path_directory_chromosomes_source=${4} # full path to source parent directory for GWAS across separate chromosomes
path_file_gwas_product=${5} # full path to product file for concatenation of GWAS summary statistics across chromosomes
path_file_frequency_product=${6} # full path to product file for concatenation of allele frequencies across chromosomes
name_directory_log_product=${7} # name for product directory for execution log files
prefix_file_log_product=${8} # file name prefix for product execution log files
suffix_file_log_product=${9} # file name suffix for product execution log files
chromosome_xy=${10} # whether to collect GWAS summary statistics, allele frequencies, and execution log reports for Chromosomes X and XY
report=${11} # whether to print reports

################################################################################
# Organize paths.

# Extract base names of files and directories.
# Removal of file suffix string using "basename" is unreliable (parameter "-s").
# "${[source string]/[string for which to replace first occurrence]/[replacement string or variable]}"
# "${[source string]//[string for which to replace all occurrences]/[replacement string or variable]}"
#name_base_file_gwas_product="$(basename $path_file_gwas_product -s .txt.gz)"
name_base_file_gwas_product="$(basename $path_file_gwas_product)"
name_base_file_gwas_product="${name_base_file_gwas_product/".txt.gz"/}"
name_base_file_frequency_product="$(basename $path_file_frequency_product)"
name_base_file_frequency_product="${name_base_file_frequency_product/".afreq.gz"/}"
path_directory_product="$(dirname $path_file_gwas_product)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_gwas_product}" # hopefully unique
path_directory_log_temporary="${path_directory_product}/${name_directory_log_product}" # hopefully unique

path_file_gwas_temporary="${path_directory_temporary}/${name_base_file_gwas_product}_concatenation.txt"
path_file_frequency_temporary="${path_directory_temporary}/${name_base_file_frequency_product}_concatenation.afreq"

# Initialize directory.
#rm -r $path_directory_product # Caution: removes the parent directory of the product files
mkdir -p $path_directory_product
rm -r $path_directory_temporary || true # silence warning if file or directory does not exist
mkdir -p $path_directory_temporary
rm -r $path_directory_log_temporary || true # silence warning if file or directory does not exist
mkdir -p $path_directory_log_temporary

# Remove any previous version of the product files.
rm $path_file_gwas_product || true # silence warning if file or directory does not exist
rm $path_file_frequency_product || true # silence warning if file or directory does not exist

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "path to source directory: " $path_directory_chromosomes_source
  echo "base file name for GWAS file: " $name_base_file_gwas_product
  echo "base file name for allele frequency file: " $name_base_file_frequency_product
  echo "path to product GWAS file: " $path_file_gwas_product
  echo "path to product frequency file: " $path_file_frequency_product
fi

################################################################################
# Execute procedure.

##########
# Initialize tables for concatenation.

# Initialize names of columns in table for concatenation of GWAS summary statistics.
# echo "#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P" > $path_gwas_concatenation
path_directory_chromosome_source="${path_directory_chromosomes_source}/chromosome_1"
#matches=("${path_source_chromosome}/${pattern_gwas_report_file}")
matches=$(find "${path_directory_chromosome_source}" -name "$pattern_file_gwas_source")
path_file_gwas_source=${matches[0]}
cat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_gwas_temporary

# Initialize names of columns in table for concatenation of allele frequencies.
path_directory_chromosome_source="${path_directory_chromosomes_source}/chromosome_1"
matches=$(find "${path_directory_chromosome_source}" -name "$pattern_file_frequency_source")
path_file_frequency_source=${matches[0]}
cat $path_file_frequency_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_frequency_temporary

##########
# Concatenation across chromosomes.

# Iterate across relevant chromosomes.
#for (( index=$chromosome_start; index<=$chromosome_end; index+=1 )); do
if [[ "$chromosome_xy" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "x" "xy")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
for chromosome in "${chromosomes[@]}"; do
  # Define path to directory for chromosome.
  path_directory_chromosome_source="$path_directory_chromosomes_source/chromosome_${chromosome}"

  # Find file for GWAS summary statistics and concatenate (collect) information.
  matches=$(find "${path_directory_chromosome_source}" -name "$pattern_file_gwas_source")
  path_file_gwas_source=${matches[0]}
  cat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 && NR < 31 { print $0 }' >> $path_file_gwas_temporary
  #cat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 { print $0 }' >> $path_file_gwas_temporary

  # Find file for allele frequencies and concatenate (collect) information.
  matches=$(find "${path_directory_chromosome_source}" -name "$pattern_file_frequency_source")
  path_file_frequency_source=${matches[0]}
  cat $path_file_frequency_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 && NR < 31 { print $0 }' >> $path_file_frequency_temporary
  #cat $path_file_frequency_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 { print $0 }' >> $path_file_frequency_temporary

  # Find and copy file for PLINK2 execution log.
  matches=$(find "${path_directory_chromosome_source}" -name "$pattern_file_log_source")
  path_file_log_source=${matches[0]}
  name_file_log_product="${prefix_file_log_product}${chromosome}${suffix_file_log_product}"
  path_file_log_product="${path_directory_log_temporary}/${name_file_log_product}"
  cp "${path_file_log_source}" "${path_file_log_product}"

done

##########
# Compress product files.

# Compress directories and files.
gzip -cvf $path_file_gwas_temporary > $path_file_gwas_product
gzip -cvf $path_file_frequency_temporary > $path_file_frequency_product
tar --remove-files -czvf "${path_directory_product}/${name_directory_log_product}.tar.gz" -C "${path_directory_log_temporary}" . # simplify path to source

##########
# Remove temporary, intermediate files.
#rm -r $path_directory_temporary
