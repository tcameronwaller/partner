#!/bin/bash

################################################################################
# Notes:

# Download newest genetic reference files:
# genomic sequences in GRCh37 and GRCh38
# dbSNP in GRCh37 and GRCh38
# Use BCFTools to annotate the reference sequences with SNP rs IDs... is that right?
# UCSC and Ensembl chain files from GRCh37 to GRCh38 and vice versa


# update installation of BCFTools and HTSlib
# install GWAS2VCF
# install BCFTools plugin for GWAS-VCF

# Steps on GWAS sum stats before SBayesR or LDPred2
# 1. convert GWAS sum stats to team standard format
# 2. convert GWAS sum stats to GWAS-VCF format
# 3. run BCFTools +Munge plugin on GWAS in GWAS-VCF format
# 4. convert from GWAS-VCF format to GWAS catalog format
# 5. convert from GWAS catalog format to format for SBayesR and LDPred2

# Steps on output from SBayesR and LDPred2
# 1. convert output to format for GWAS2VCF
# 2. convert to GWAS-VCF format
# 3. run BCFTools +Liftover plugin from GRCh37 to GRCh38
# 4. convert from GWAS-VCF format to GWAS catalog format
# 5. convert from GWAS catalog format to format readable by PLINK2 score function


################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
path_environment_gwas2vcf="${path_directory_tools}/python/environments/gwas2vcf"
path_gwas2vcf="${path_directory_tools}/gwas2vcf/gwas2vcf/main.py"

#path_plink2=$(<"./tools_plink2.txt")
#path_gctb=$(<"./tools_waller_gctb.txt")
#path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
#path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="${path_directory_dock}/test_gwas_clean"
path_directory_reference_gwas2vcf="${path_directory_dock}/test_gwas_clean/reference_gwas2vcf"

# Files.
path_file_gwas_standard_source="${path_directory_dock}/hormone_genetics/gwas_format_standard/32042192_ruth_2020_testosterone_female.txt.gz"
path_file_gwas_product="${path_directory_product}/32042192_ruth_2020_testosterone_female.txt.gz"
path_file_gwas2vcf_parameter="${path_directory_product}/importation_parameter_gwas_standard_to_gwasvcf.json"

# Scripts.
#path_script_gwas_format_source_to_standard="${path_directory_process}/promiscuity/scripts/gwas_format/format_gwas_team/translate_gwas_30367059_teumer_2018.sh"
path_script_access_reference_gwas2vcf="${path_directory_process}/promiscuity/scripts/gwas_clean/access_reference_gwas2vcf.sh"
#path_script_gwas_format_standard_to_gwas2vcf="${path_directory_process}/promiscuity/scripts/utility/gwas_clean/translate_gwas_standard_to_gwas2vcf.sh" # <-- do not need
#path_script_gwas_clean="${path_directory_process}/promiscuity/scripts/utility/gwas_clean/clean_gwas.sh"

# Initialize directories.
#rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



################################################################################
# Organize parameters.


report="true"

################################################################################
# Execute procedure.



##########
# Translate GWAS summary statistics to standard format.
# Product Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT
if false; then
  /usr/bin/bash "${path_script_gwas_format}" \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi



##########
# Installation: GWAS2VCF
# PubMed: 33441155
# GWAS-VCF format specification: https://github.com/MRCIEU/gwas-vcf-specification
# Host of GWAS2VCF: https://github.com/MRCIEU/gwas2vcf
# Documentation for GWAS2VCF: https://mrcieu.github.io/gwas2vcf
# Refer to notes in script "install_local_software.sh".
# Refer to notes in script "install_python_virtual_environments_packages.sh"
if false; then
  # Activate Virtual Environment.
  source "${path_environment_gwas2vcf}/bin/activate"
  echo "confirm Python Virtual Environment path..."
  which python3
  sleep 5s
  # Test installation of GWAS2VCF.
  python3 $path_gwas2vcf -h
  # Deactivate Virtual Environment.
  deactivate
  which python3
fi



##########
# Access genomic reference information for GWAS2VCF.

if true; then
  /usr/bin/bash $path_script_access_reference_gwas2vcf \
  $path_directory_reference_gwas2vcf \
  $report
fi



##########
# Translate GWAS summary statistics from standard format to GWAS-VCF format.
# Use tool GWAS2VCF.
# Functions.
# 1. Verify and introduce SNPs' rs identifiers from dbSNP reference.


# Indexing of columns in source GWAS summary statistics bases on zero.
# TODO: decompress the summary statistics
if false; then
  echo "{
    "chr_col": 1,
    "pos_col": 2,
    "snp_col": 0,
    "ea_col": 3,
    "oa_col": 4,
    "beta_col": 6,
    "se_col": 7,
    "ncontrol_col": 9,
    "pval_col": 8,
    "eaf_col": 5,
    "delimiter": "\t",
    "header": true,
    "build": "GRCh37"
  }" > $path_file_gwas2vcf_parameter
fi




##########
# Munge GWAS summary statistics in Bioconductor package "MungeSumstats".
# Documentation: https://bioconductor.org/packages/release/bioc/manuals/MungeSumstats/man/MungeSumstats.pdf



##########
# TCW; 24 January 2023
# The Bioconductor packange "MungeSumstats" might be preferrable over the
# "+munge" plugin for BCFTools.
##########
# Munge GWAS summary statistics in BCFTools "+munge" plugin.
# documentation: https://github.com/freeseek/score#convert-summary-statistics


#
