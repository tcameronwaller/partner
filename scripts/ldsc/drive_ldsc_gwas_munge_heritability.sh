#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...
# ...
# TODO: TCW; 14 July 2022
# TODO: 1. separate script for LDSC Munge
# TODO: 2. separate script for LDSC heritability (from GWAS sum-stats already munged)
# TODO: 3. separate script for LDSC correlation (from GWAS sum-stats already munged)
# TODO: 4. LDSC procedures use NumPy which seems to use as many threads as possible
# TODO: - - restrict threads in all LDSC driver scripts...
# TODO: - - use "MKL_NUM_THREADS" global variable to limit threads
# TODO: https://groups.google.com/g/ldsc_users/c/Bnj7FFl5jlw/m/nrQ7yH-8BgAJ?pli=1

# Force SciPy not to use all available cores on a cluster computation node.
#export MKL_NUM_THREADS=$threads

################################################################################
################################################################################
################################################################################

# TODO: Need to include liability scale in heritability estimates for case-control studies

################################################################################
# Organize arguments.
path_gwas_source_parent=${1} # full path to parent directory for source GWAS summary statistics
path_gwas_target_parent=${2} # full path to parent directory for target GWAS summary statistics
path_heritability_parent=${3} # full path to directory for heritability report
path_genetic_reference=${4} # full path to directory for genetic reference information
response=${5} # whether GWAS response is beta coefficient ("coefficient"), odds ratio ("odds_ratio"), or z-scores ("z_score")
report=${6} # whether to print reports

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
# Argument "--a1-inc" tells LDSC that GWAS sumstats do not have signed statistics but are coded so that A1 allele always increases effect.
# Argument ""--signed-sumstats" designates the column name for signed statistics.
# Examples:
# - "--signed-sumstats BETA,0"
# - "--signed-sumstats OR,1"
# - "--signed-sumstats Z,0"

if [[ "$response" == "coefficient" ]]; then
  $path_ldsc/munge_sumstats.py \
  --sumstats $path_gwas_format_compress \
  --signed-sumstats BETA,0 \
  --merge-alleles $path_alleles/w_hm3.snplist \
  --out $path_gwas_munge
elif [[ "$response" == "z_score" ]]; then
  $path_ldsc/munge_sumstats.py \
  --sumstats $path_gwas_format_compress \
  --signed-sumstats Z,0 \
  --merge-alleles $path_alleles/w_hm3.snplist \
  --out $path_gwas_munge
elif [[ "$response" == "odds_ratio" ]]; then
  $path_ldsc/munge_sumstats.py \
  --sumstats $path_gwas_format_compress \
  --signed-sumstats OR,1 \
  --merge-alleles $path_alleles/w_hm3.snplist \
  --out $path_gwas_munge
else
  echo "invalid specification of GWAS effect"
fi

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
