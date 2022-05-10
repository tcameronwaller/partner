#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...
# 1. I still need to use "git clone" to pull down latest release versions
# of PRS-CS and PRS-CSX to the server
# 2. I then need to specify paths to those local source copies of PRS-CS and PRS-CSX
# 3. I also need some sort of "validation .bim" file for a list of relevant SNPs
# - - these SNPs must come from the 'target' genotypes for the Polygenic Scores
# 4. I also need to format the PLINK2 linear and logistic GWAS reports for PRS-CS format.
################################################################################
################################################################################
################################################################################


################################################################################
# Organize arguments.
path_genetic_reference=${1} # full path to directory for genetic references
path_source_gwas_summary=${2} # full path to source GWAS summary statistics in format for PRS-CS
count_gwas_samples=${3} # integer count of samples for the GWAS study
population_ancestry=${4} # character code of ancestral population of GWAS study
path_target_directory=${5} # full path to directory for product reports
file_name_prefix=${6} # prefix for names of product report files
report=${7} # whether to print reports

###########################################################################
# Organize paths.

path_1000_genomes="$path_genetic_reference/1000_genomes"
path_uk_biobank="$path_genetic_reference/uk_biobank"

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
path_prs_cs=$(<"./tools_prs_cs.txt")
path_prs_csx=$(<"./tools_prs_csx.txt")
path_environment_prs_cs="${path_tools}/python/environments/prs_cs"
source "${path_environment_prs_cs}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s

################################################################################
# Test PRS-CSX.
# https://github.com/getian107/PRScsx

python3 $path_prs_csx/PRScsx.py \
--ref_dir=$path_1000_genomes \
--bim_prefix=$path_bim_snp_relevance \
--sst_file=$path_source_gwas_summary \
--n_gwas=$count_gwas_samples \
--pop=$population_ancestry \
--out_dir=$path_target_directory \
--out_name=$file_name_prefix \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_heritability_report

################################################################################
# Deactivate Virtual Environment.
deactivate
which python3

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "drive_prscs_test.sh"
  echo "----------"
  echo "PRS-CS report:"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
