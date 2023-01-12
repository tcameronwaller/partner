#!/bin/bash

################################################################################
################################################################################
################################################################################
# Notes...

################################################################################
################################################################################
################################################################################

################################################################################
# Organize arguments.

path_file_gwas=${1} # full directory path and file name for source GWAS summary statistics in GCTB and GCTA-COJO ".ma" format without compression
path_file_ld_matrix=${2} # full directory path and base file name for Linkage Disequilibrium (LD) reference matrix in GCTB format
path_file_base_product=${3} # full directory path and base file name for product files from GCTB SBayesR
path_gctb=${4} # full directory path and file name for local executable installation of GCTB SBayesR
report=${5} # whether to print reports

################################################################################
# Run SBayesR.

##########
# Apply SBayesR to adjust weights of effect sizes across SNPs.
# The path to the LD matrix actually points to three separate files with
# different suffixes: ".bin", ".info", ".log".
# Extra commands for SBayesR:
# --unscale-genotype
# --exclude-mhc
# --exclude-region
# --impute-n

$path_gctb \
--sbayes R \
--exclude-mhc \
--ldm $path_file_ld_matrix \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--gwas-summary $path_file_gwas \
--chain-length 10000 \
--burn-in 2000 \
--out-freq 10 \
--out $path_file_base_product



#
