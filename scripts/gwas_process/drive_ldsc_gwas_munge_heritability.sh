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
path_gwas_source_parent=${1} # full path to parent directory for source GWAS summary statistics
path_gwas_target_parent=${2} # full path to parent directory for target GWAS summary statistics
path_heritability_parent=${3} # full path to directory for heritability report
path_genetic_reference=${4} # full path to directory for genetic reference information
report=${4} # whether to print reports

################################################################################
# Paths.
path_alleles="${path_genetic_reference}/alleles"
path_disequilibrium="${path_genetic_reference}/disequilibrium"
path_baseline="${path_genetic_reference}/baseline"
path_weights="${path_genetic_reference}/weights"
path_frequencies="${path_genetic_reference}/frequencies"

path_gwas_format_compress="${path_gwas_source_parent}/gwas_format.txt.gz"
path_gwas_munge="${path_gwas_target_parent}/gwas_munge"
path_gwas_munge_suffix="${path_gwas_target_parent}/gwas_munge.sumstats.gz"
path_gwas_munge_log="${path_gwas_target_parent}/gwas_munge.log"
path_heritability_report="${path_heritability_parent}/heritability_report"
path_heritability_report_suffix="${path_heritability_parent}/heritability_report.log"

# Remove any previous versions of target files.
rm $path_gwas_munge_suffix
rm $path_gwas_munge_log
rm $path_heritability_report_suffix

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

################################################################################
# Munge GWAS summary statistics for use in LDSC.
$path_ldsc/munge_sumstats.py \
--sumstats $path_gwas_format_compress \
--out $path_gwas_munge \
--merge-alleles $path_alleles/w_hm3.snplist
#--a1-inc # This flag tells LDSC that GWAS sumstats do not have signed Betas but are coded so that A1 allele always increases effect
#--signed-sumstats BETA,0 # Brandon J. Coombes used this argument in a script, but I don't think it is necessary

################################################################################
# Estimate phenotype heritability in LDSC.
$path_ldsc/ldsc.py \
--h2 $path_gwas_munge_suffix \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_heritability_report

################################################################################
# Deactivate Virtual Environment.
deactivate
which python2

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_ldsc_gwas_munge_heritability.sh"
  echo "----------"
  echo "LDSC heritability report:"
  cat $path_heritability_report_suffix
  echo "----------"
fi
