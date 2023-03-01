#!/bin/bash

################################################################################
# Notes:

# This script accesses Linkage Disequilibrium (LD) reference matrices for use in
# GCTB and SBayesR.

# GCTB host website.
# Site: https://cnsgenomics.com/software/gctb/#Download

# Review: TCW; 1 March 2023
# Last Accession: TCW; 1 March 2023 (to NCSA mForge server)

################################################################################



################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

# Initialize directory.
mkdir -p $path_directory_parent

cd $path_directory_parent

###########################################################################
# Execute procedure.

# Echo each command to console.
#set -x
# Suppress echo each command to console.
#set +x

###########################################################################
# Organize directories.
# Access reference information from NCBI dbSNP.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Accessing Linkage Disequilibrium (LD) reference matrices for use in"
  echo "GCTB SBayesR."
  echo "Human Genome Assembly: GRCh37"
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Access the specific files and save within the directory.

# GCTB SBayesR "Shrunk sparse matrix"
# PubMed: 31704910
# Human Genome Assembly: GRCh37
# Variants: 1.1 million HapMap3 SNPs
# Cohort: 50,000 randomly-select unrelated persons of European ancestry from UK Biobank
# Detail: observed LD correlations shrunk toward expected value from genetic map on 1000 Genomes (PubMed: 21479081); LD correlations less than threshold (1e-5) set to zero for sparse format;
# Host site: https://zenodo.org/record/3350914
# File date: 23 August 2019; File size: 22.1 GB
if true; then
  wget --directory-prefix $path_directory_parent --content-disposition --no-check-certificate "https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip"
fi

# GCTB SBayesR "Sparse matrix (including MHC regions)"
# PubMed: 33608517
# Human Genome Assembly: GRCh37
# Variants: 1.1 million HapMap3 SNPs
# Cohort: 50,000 unrelated persons of European ancestry from UK Biobank
# Detail: LD correlations set to zero according to chi-squared test statistic threshold (10)
# Host site: https://cnsgenomics.com/software/gctb/#Download
# File date: ?; File size: ___ GB
if false; then
  wget --directory-prefix $path_directory_parent --content-disposition --no-check-certificate "https://cnsgenomics.com/data/GCTB/ukbEURu_imp_v3_HM3_n50k.chisq10.zip"
fi

# GCTB SBayesR "Banded matrix (including MHC regions)"
# PubMed: 33326037
# Human Genome Assembly: GRCh37
# Variants: 1.1 million HapMap3 SNPs
# Cohort: 10,000 unrelated persons of European ancestry from UK Biobank
# Detail: LD correlation banded matrix with window size of 3 cM per SNP (PubMed: 33326037);
# Host site: https://cnsgenomics.com/software/gctb/#Download
# File date: ?; File size: __ GB
if false; then
  wget --directory-prefix $path_directory_parent --content-disposition --no-check-certificate "https://cnsgenomics.com/data/GCTB/band_ukb_10k_hm3.zip"
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "access_gctb_reference_ld_matrices.sh"
  echo "----------"
fi





#
