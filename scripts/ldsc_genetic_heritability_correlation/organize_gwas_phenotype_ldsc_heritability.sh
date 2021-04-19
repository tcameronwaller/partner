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
study=${1} # identifier of GWAS study
source_file=${2} # name of source file with GWAS summary statistics
path_source_file=${3} # full path to source file with GWAS summary statistics
path_genetic_reference=${4} # full path to parent directory with genetic reference files for LDSC
path_gwas=${5} # full path to parent directory for formatted GWAS summary statistics
path_heritability=${6} # full path to parent directory for LDSC heritability estimation
path_genetic_correlation=${7} # full path to parent directory for LDSC genetic correlation
path_script_gwas_format=${8} # full path to script to use to organize format of GWAS summary statistics for phenotype
path_promiscuity_scripts=${9} # complete path to directory of scripts for z-score standardization
report=${10} # whether to print reports

################################################################################
# Organize variables.

# Read private, local file paths.
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")

path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"

path_study_gwas="${path_gwas}/${study}"
path_study_heritability="${path_heritability}/${study}"
path_study_genetic_correlation="${path_genetic_correlation}/${study}"
mkdir -p $path_study_gwas
mkdir -p $path_study_heritability
mkdir -p $path_study_genetic_correlation
path_gwas_collection="${path_study_gwas}/gwas_collection.txt"
path_gwas_format="${path_study_gwas}/gwas_format.txt"
path_gwas_format_compress="${path_gwas_format}.gz"
path_gwas_munge="${path_study_gwas}/gwas_munge"
path_gwas_munge_suffix="${path_gwas_munge}.sumstats.gz"
path_gwas_munge_log="${path_gwas_munge}.log"
path_heritability_report="${path_study_heritability}/heritability_report"
path_heritability_report_suffix="${path_heritability_report}.log"

#path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_4_of_5.sh"
path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_5_of_6.sh"

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------"
  echo "file name: " $source_file
  echo "path to original file: " $path_source_file
  echo "path to new file: " $path_gwas_format
  echo "----------"
fi

# Organize information in format for LDSC.
# Parameters.
/usr/bin/bash "$path_script_gwas_format" \
$study \
$path_source_file \
$path_gwas_collection \
$path_gwas_format \
$path_gwas_format_compress \
$path_promiscuity_scripts \
$report

# Munge GWAS summary statistics for use in LDSC.
$path_ldsc/munge_sumstats.py \
--sumstats $path_gwas_format_compress \
--out $path_gwas_munge \
--merge-alleles $path_alleles/w_hm3.snplist \
#--a1-inc

# Heritability.
$path_ldsc/ldsc.py \
--h2 $path_gwas_munge_suffix \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_heritability_report

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "LDSC heritability report:"
  cat $path_heritability_report_suffix
  echo "----------"
fi

###########################################################################
# Remove temporary files.
rm $path_gwas_collection
rm $path_gwas_format
