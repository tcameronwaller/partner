#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

# TODO: re-write this... define the study_gwas and study_heritability directories upstream...
# That will make it possible to use this for metabolites...


################################################################################
# Organize arguments.
study=${1} # identifier of GWAS study
name_prefix=${2} # file name prefix for GWAS files and heritability file or "null"
path_source_file=${3} # full path to source file with GWAS summary statistics
path_genetic_reference=${4} # full path to parent directory with genetic reference files for LDSC
path_gwas_study=${5} # full path to parent directory for formatted GWAS summary statistics
path_heritability_study=${6} # full path to parent directory for LDSC heritability estimation
path_script_gwas_format=${7} # full path to script to use to organize format of GWAS summary statistics for phenotype
path_promiscuity_scripts=${8} # complete path to directory of scripts for z-score standardization
path_ldsc=${9} # path to LDSC
report=${10} # whether to print reports

################################################################################
# Organize variables.

path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"

if [[ "$name_prefix" != "null" ]]; then
  path_gwas_collection="${path_gwas_study}/${name_prefix}_gwas_collection.txt"
  path_gwas_format="${path_gwas_study}/${name_prefix}_gwas_format.txt"
  path_gwas_munge="${path_gwas_study}/${name_prefix}_gwas_munge"
  path_heritability_report="${path_heritability_study}/${name_prefix}_heritability_report"
fi

if [[ "$name_prefix" == "null" ]]; then
  path_gwas_collection="${path_gwas_study}/gwas_collection.txt"
  path_gwas_format="${path_gwas_study}/gwas_format.txt"
  path_gwas_munge="${path_gwas_study}/gwas_munge"
  path_heritability_report="${path_heritability_study}/heritability_report"
fi

path_gwas_format_compress="${path_gwas_format}.gz"
path_gwas_munge_suffix="${path_gwas_munge}.sumstats.gz"
path_gwas_munge_log="${path_gwas_munge}.log"
path_heritability_report_suffix="${path_heritability_report}.log"

#path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_4_of_5.sh"
path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_5_of_6.sh"

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------"
  echo "path to original source file: " $path_source_file
  echo "path to new file: " $path_gwas_format
  echo "----------"
fi

# Initialize directories.
if [ ! -d $path_gwas_study ]; then
    mkdir -p $path_gwas_study
fi
if [ ! -d $path_heritability_study ]; then
    mkdir -p $path_heritability_study
fi

# Organize information in format for LDSC.
# Parameters.
rm $path_gwas_format
rm $path_gwas_format_compress
/usr/bin/bash "$path_script_gwas_format" \
$study \
$path_source_file \
$path_gwas_collection \
$path_gwas_format \
$path_gwas_format_compress \
$path_promiscuity_scripts \
$report

# Munge GWAS summary statistics for use in LDSC.
rm $path_gwas_munge_suffix
rm $path_gwas_munge_log
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
#rm $path_gwas_munge_suffix
#rm $path_gwas_munge_log
