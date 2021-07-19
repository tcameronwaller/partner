#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "Munge GWAS summary statistics for LDSC."
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo ""
echo ""
echo ""

# Organize variables.
phenotype_one=${1} # name of phenotype one
phenotype_two=${2} # name of phenotype two
path_gwas_one=${3} # path to GWAS report files for phenotype one
path_gwas_two=${4} # path to GWAS report files for phenotype two
path_gwas_munge_one=${5} # path for munge GWAS report files for phenotype one
path_gwas_munge_two=${6} # path for munge GWAS report files for phenotype two
path_report=${7} # path to directory for genetic correlation report
path_alleles=${8} # path to reference for alleles
path_disequilibrium=${9} # path to reference for linkage disequilibrium
path_ldsc=${10} # path to LDSC

###########################################################################
# Munge GWAS summary statistics for LDSC.

# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
# "(Note that some summary statistic files do not have a signed summary
# statistic, but are coded so that A1 is always the trait- or risk-increasing
# allele. This is equivalent to providing a signed summary statistic, and
# munge_sumstats.py will process such files if called with the --a1-inc1 flag)."

# PLINK2 reports Odds Ratios (logistic regression) and Beta Coefficients
# (linear regression) relative to the A1 allele.

# I think that it is not necessary to use the "--a1-inc" flag.

path_gwas_munge_one_suffix="${path_gwas_munge_one}.sumstats.gz"
path_gwas_munge_one_log="${path_gwas_munge_one}.log"
$path_ldsc/munge_sumstats.py \
--sumstats $path_gwas_one \
--out $path_gwas_munge_one \
--merge-alleles $path_alleles/w_hm3.snplist \
#--a1-inc

path_gwas_munge_two_suffix="${path_gwas_munge_two}.sumstats.gz"
path_gwas_munge_two_log="${path_gwas_munge_two}.log"
$path_ldsc/munge_sumstats.py \
--sumstats $path_gwas_two \
--out $path_gwas_munge_two \
--merge-alleles $path_alleles/w_hm3.snplist \
#--a1-inc

###########################################################################
# Estimate genetic correlation in LDSC.

$path_ldsc/ldsc.py \
--rg $path_gwas_munge_one_suffix,$path_gwas_munge_two_suffix \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_report
