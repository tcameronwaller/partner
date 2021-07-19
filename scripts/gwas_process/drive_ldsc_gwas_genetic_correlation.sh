#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...
# ...
################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.
path_gwas_primary_munge_suffix=${1} # full path to parent directory for source GWAS summary statistics
path_gwas_secondary_munge_suffix=${2} # full path to parent directory for target GWAS summary statistics
path_genetic_correlation_parent=${3} # full path to directory for genetic correlation report
path_genetic_reference=${4} # full path to directory for genetic reference information
report=${5} # whether to print reports

################################################################################
# Paths.
path_alleles="${path_genetic_reference}/alleles"
path_disequilibrium="${path_genetic_reference}/disequilibrium"
path_baseline="${path_genetic_reference}/baseline"
path_weights="${path_genetic_reference}/weights"
path_frequencies="${path_genetic_reference}/frequencies"

path_correlation_report="${path_genetic_correlation_parent}/correlation_report"
path_correlation_report_suffix="${path_genetic_correlation_parent}/correlation_report.log"

# Remove any previous versions of target files.
rm $path_correlation_report_suffix

################################################################################
# Activate Virtual Environment.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_ldsc=$(<"./tools_ldsc.txt")
path_environment_ldsc="${path_tools}/python/environments/ldsc"
source "${path_environment_ldsc}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python2
sleep 5s

echo "here is the out path"
echo $path_genetic_correlation_report

################################################################################
# Estimate genetic correlation in LDSC.
$path_ldsc/ldsc.py \
--rg ${path_gwas_primary_munge_suffix},${path_gwas_secondary_munge_suffix} \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_genetic_correlation_report

################################################################################
# Deactivate Virtual Environment.
deactivate
which python2

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_ldsc_gwas_genetic_correlation.sh"
  echo "----------"
  echo "LDSC genetic correlation report:"
  cat $path_correlation_report_suffix
  echo "----------"
fi
