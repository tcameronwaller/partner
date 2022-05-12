#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...

# This script calls PRS-CS to generate posterior effects across SNP alleles.

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.
path_source_gwas_summary=${1} # full path to source GWAS summary statistics in format for PRS-CS
count_gwas_samples=${2} # integer count of samples for the GWAS study
path_snp_relevance_bim_container=${3} # full path to directory for files with relevant SNPs in '.bim' format
path_target_directory_score=${4} # full path to directory for product reports on posterior effect size estimates
file_name_prefix=${5} # prefix for names of product report files
path_genetic_reference_prscs=${6} # full path to directory for genetic references
population_ancestry=${7} # character code of ancestral population of GWAS study: 'AFR', 'AMR', 'EAS', 'EUR', or 'SAS'
chromosome_x=${8} # whether to collect GWAS summary statistics report for Chromosome X
path_promiscuity_scripts=${9} # full path to directory of general scripts
report=${10} # whether to print reports

###########################################################################
# Organize paths.

path_1000_genomes="${path_genetic_reference_prscs}/1000_genomes"
path_uk_biobank="${path_genetic_reference_prscs}/uk_biobank"

# Scripts.
path_scripts_gwas_process="${path_promiscuity_scripts}/gwas_process"
path_script_estimate_prscs="${path_scripts_gwas_process}/prscs_polygenic_score/estimate_prscs_allelic_effects.sh"

# Determine relevant chromosomes.
if [[ "$chromosome_x" == "true" ]]; then
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
else
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
fi
for chromosome in "${chromosomes[@]}"; do
  # Define paths for chromosome.
  name_bim_file_chromosome="MERGED.maf0.dosR20.3.noDups.chr${chromosome}.dose.vcf.gz" # omit the '.bim' suffix
  path_snp_relevance_bim_prefix="$path_snp_relevance_bim_container/${name_bim_file_chromosome}"
  # Define name of product file for chromosome.
  name_file_product="${file_name_prefix}_chromosome${chromosome}"
  # Call script.
  /usr/bin/bash "${path_script_estimate_prscs}" \
  $path_source_gwas_summary \
  $count_gwas_samples \
  $path_snp_relevance_bim_prefix \
  $path_target_directory_score \
  $name_file_product \
  $path_1000_genomes \
  $population_ancestry \
  $chromosome \
  $report
done
