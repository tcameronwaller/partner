#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...

# This script WILL CALL PLINK2 to apply a linear combination across allelic
# effects to create a polygenic score to represent an entire genotype.

################################################################################
################################################################################
################################################################################

# 1. I need a "$path_genotype" file


################################################################################
# Organize arguments.
path_source_prscs_allelic_effects=${1} # full path to PRS-CS allelic effects with concatenation across chromosomes
report=${2} # whether to print reports

###########################################################################
# Organize paths.

################################################################################
# Activate Virtual Environment.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_plink2=$(<"./tools_plink2.txt")

################################################################################
# Combine allelic effects in PLINK2.
# https://www.cog-genomics.org/plink/1.9/score
# https://www.cog-genomics.org/plink/2.0/score

# Format from PRS-CS.
#   chromosome   rsID   base position   A1   A2   effect size


$path_plink2 \
--memory 90000 \
--threads $threads \
--bgen $path_genotype ref-first \
--sample $path_sample \
--xchr-model 2 \
--score $path_source_prscs_allelic_effects 2 4 header-read no-mean-imputation ignore-dup-ids \
--score-col-nums 6 \
--out report

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_prscs_plink_polygenic_linear_combination_score.sh"
  echo "----------"
  echo "___ report:"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
