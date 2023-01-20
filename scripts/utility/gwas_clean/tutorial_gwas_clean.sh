#!/bin/bash

################################################################################
# Notes:

# Download newest genetic reference files:
# genomic sequences in GRCh37 and GRCh38
# dbSNP in GRCh37 and GRCh38
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
