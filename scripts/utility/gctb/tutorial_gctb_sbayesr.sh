#!/bin/bash

################################################################################
# Notes:
# This script includes notes and code to work through an introductory tutorial
# on using the GCTB tool package.
# GCTB is a package of tools for Genome-wide Complex Trait Bayesian (GCTB)
# analysis that includes the tool SBayesR for calculation of Polygenic Scores
# (PGS).
# Documentation on GCTB.
# site: https://cnsgenomics.com/software/gctb/#Overview
# tutorial: https://cnsgenomics.com/software/gctb/#Tutorial
# educational presentation slides and video: https://www.cnsgenomics.com/data/teaching/PCTG/SBayes/

# Author: T. Cameron Waller
# Date, first execution: 10 January 2023
# Date, review: ____

# It is important that Genome-wide Association Study (GWAS) summary statistics
# and Linkage Disequilibrium (LD) reference matrices to use names (rs
# identifiers) and coordinates for Single Nucleotide Polymorphisms (SNPs) that
# are in the same assembly of the human genome, such as GRCh37 or GRCh38.
# As of year 2023, LD reference matrices for GCTB SBayesR used GRCh37.

# Note: TCW; 13 January 2023
# SBayesR does not converge when using the tutorial LD matrix (chromosome 22
# only) with real GWAS summary statistics, even after filtering the summary
# statistics to chromosome 22 only. A possible explanation would be that the
# coefficients (betas) in the GWAS summary statistics would not have the same
# distribution after filtering to chromosome 22 only.

# Note: TCW; 17 January 2023
# In a test on real GWAS summary statistics and the shrunk sparse LD matrix for
# chromosome 1 (largest chromosome), the SBayesR procedure completed
# successfully in fewer than 10 minutes on the NCSA mForge head node 2.

# Note: TCW; 17 January 2023
# A test running GCTB SBayesR procedure on all autosomal chromosomes (1-22) for
# a set of real GWAS summary statistics completed successfully.

# TODO: TCW; 17 January 2023
# Now I need to translate SNP names (rsID), and genomic coordinates to GRCh38.


################################################################################

##########
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
path_plink2=$(<"./tools_plink2.txt")
path_gctb=$(<"./tools_waller_gctb.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parent="${path_directory_dock}/test_sbayesr"

# Files.
path_file_gwas_source="${path_directory_dock}/hormone_genetics/gwas_format_standard/32042192_ruth_2020_testosterone_bioavailable_female.txt.gz"
path_file_gwas_source_temporary="${path_directory_parent}/32042192_ruth_2020_testosterone_bioavailable_female_temp.txt"
path_file_gwas_source_temporary_compress="${path_directory_parent}/32042192_ruth_2020_testosterone_bioavailable_female_temp.txt.gz"
path_file_gwas_product="${path_directory_parent}/32042192_ruth_2020_testosterone_bioavailable_female.ma"
path_file_gwas_tutorial="${path_directory_parent}/gctb_2.0_tutorial/ma/sim_1.ma"
path_file_ld_matrix_tutorial="${path_directory_parent}/gctb_2.0_tutorial/ldm/sparse/chr22/1000G_eur_chr22.ldm.sparse"
path_file_ld_matrix_chromosome_1="${path_directory_parent}/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr1_v3_50k.ldm.sparse"
path_file_base_product="${path_directory_parent}/test_sbayesr_out_tcw_723"

# Scripts.
path_script_gwas_format="${path_directory_process}/promiscuity/scripts/utility/gctb/constrain_translate_gwas_standard_to_gctb.sh"
path_script_run_sbayesr="${path_directory_process}/promiscuity/scripts/utility/gctb/run_gctb_sbayesr.sh"

# Uniform Resource Locators (URLs).
url_gctb="https://cnsgenomics.com/software/gctb/download/gctb_2.04.3_Linux.zip"
url_gctb_tutorial="https://cnsgenomics.com/software/gctb/download/gctb_2.0_tutorial.zip"
# Linkage Disequilibrium (LD) matrices for reference:
# 1. Sparse, shrinkage LD matrix on 50K people from UK Biobank used in Lloyd-Jones et al, Nature Communications, 2019 (PubMed:31704910).
url_ld_matrix_sparse_shrinkage_hapmap3="https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip"
# 2. Sparse, shrinkage LD matrix used in Zeng et al, Nature Communications, 2021 (PubMed:33608517)
url_ld_matrix_sparse_chisquare="https://cnsgenomics.com/data/GCTB/ukbEURu_imp_v3_HM3_n50k.chisq10.zip"
# 3. Banded matrix, window size 3 cM per SNP as described in Prive et al, Bioinformatics, 2020 (PubMed:33326037)
url_ld_matrix_banded="https://cnsgenomics.com/data/GCTB/band_ukb_10k_hm3.zip"



# Initialize directories.
###rm -r $path_directory_parent # Removing parent directory would lose genetic reference data.
mkdir -p $path_directory_parent
cd $path_directory_parent



###########################################################################
# Organize parameters.

report="true"



###########################################################################
# Procedure.



##########
# Installation.
# Authors provided a "statically linked 64-bit Linux executable, gctb".
# It is practical to run GCTB by a simple execution call to this handle.
# The authors also gave more thorough instructions for installation from source
# with local custom, local compilation.
# site: https://cnsgenomics.com/software/gctb/download/README.html
#cd ./ # navigate to directory in which to install GCTB.
if false; then
  # Access the specific file and save within the child directory.
  wget \
  "${url_gctb}" \
  --directory-prefix "${path_directory_tools}" \
  --content-disposition \
  --show-progress
  # Decompress content from Zip archive format.
  unzip "${path_directory_tools}/gctb_2.04.3_Linux.zip" -d "${path_directory_tools}"
  "${path_directory_tools}/gctb_2.04.3_Linux/gctb"
  $path_sbayesr
fi



##########
# Accession of data.

# Data for tutorial.
if false; then
  # Access the specific file and save within the child directory.
  wget \
  "${url_gctb_tutorial}" \
  --directory-prefix "${path_directory_parent}" \
  --content-disposition \
  --show-progress
  # Decompress content from Zip archive format.
  unzip "${path_directory_parent}/gctb_2.0_tutorial.zip" -d "${path_directory_parent}"
fi

# Linkage Disequilibrium (LD) reference matrix.
# Sparse (shrinkage) LD reference matrix on 50,000 persons of European ancestral
# background from the UK Biobank.
if false; then
  # Access the specific file and save within the child directory.
  wget \
  "${url_ld_matrix_sparse_shrinkage_hapmap3}" \
  --directory-prefix "${path_directory_parent}" \
  --content-disposition \
  --show-progress
  # Decompress content from Zip archive format.
  unzip "${path_directory_parent}/ukbEURu_hm3_sparse.zip" -d "${path_directory_parent}"
fi



##########
# Prepare GWAS Summary Statistics.

# Filter GWAS summary statistics to Chromosome 22 and translate to format for GCTB.
# For test using the tutorial LD Matrix for Chromosome 22.
if false; then
  zcat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_gwas_source_temporary
  zcat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
    if ( NF != 14)
      # Skip any rows with incorrect count of column fields.
      next
    else if ( ( $2 != "NA" ) && ( ($2 + 0) == 22 ) )
      # Print the row entirely.
      print $0
    else
      # Skip the row.
      next
    }' >> $path_file_gwas_source_temporary

  # Compress file format.
  gzip -cvf $path_file_gwas_source_temporary > $path_file_gwas_source_temporary_compress

  # Translate GWAS summary statistics.
  /usr/bin/bash "${path_script_gwas_format}" \
  $path_file_gwas_source_temporary_compress \
  $path_file_gwas_product \
  $report

  # Remove temporary files.
  rm $path_file_gwas_source_temporary
  rm $path_file_gwas_source_temporary_compress
fi

# Translate GWAS summary statistics to format for GCTB.
# Translate GWAS summary statistics to format for GCTB.
# columns: SNP   A1   A2   freq   b   se   p   N
if false; then
  /usr/bin/bash "${path_script_gwas_format}" \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi


##########
# Apply SBayesR to adjust weights of effect sizes across SNPs.
# The path to the LD matrix actually points to three separate files with
# different suffixes: ".bin", ".info", ".log".
# Extra commands for SBayesR:
# --unscale-genotype
# --exclude-mhc
# --exclude-region
# --impute-n
# File with suffix ".snpRes" gives the new effect sizes across SNPs after
# adjustment of weights for LD (I think; TCW; 12 January 2023).

# Test using tutorial GWAS summary statistics and tutorial LD matrix.
if false; then
  $path_gctb \
  --sbayes R \
  --exclude-mhc \
  --ldm $path_file_ld_matrix_tutorial \
  --pi 0.95,0.02,0.02,0.01 \
  --gamma 0.0,0.01,0.1,1 \
  --gwas-summary $path_file_gwas_tutorial \
  --chain-length 10000 \
  --burn-in 2000 \
  --out-freq 10 \
  --out $path_file_base_product 2>&1 | tee "${path_file_base_product}.log"
fi

# Test using real GWAS summary statistics and LD matrix for chromosome 1.
if true; then
  $path_gctb \
  --sbayes R \
  --exclude-mhc \
  --ldm $path_file_ld_matrix_chromosome_1 \
  --pi 0.95,0.02,0.02,0.01 \
  --gamma 0.0,0.01,0.1,1 \
  --gwas-summary $path_file_gwas_product \
  --chain-length 10000 \
  --burn-in 2000 \
  --out-freq 10 \
  --out $path_file_base_product
fi


if false; then
  /usr/bin/bash $path_script_run_sbayesr \
  $path_file_gwas_tutorial \
  $path_file_ld_matrix_tutorial \
  $path_file_base_product \
  $path_gctb \
  $report
fi

if false; then
  /usr/bin/bash $path_script_run_sbayesr \
  $path_file_gwas_product \
  $path_file_ld_matrix_tutorial \
  $path_file_base_product \
  $path_gctb \
  $report
fi




##########
# Translate SNP effects from GRCh37 to GRCh38 for calculation of polygenic
# scores on genotypes.


##########
# Apply PLINK to calculate scores across SNPs in genotypes.






#
