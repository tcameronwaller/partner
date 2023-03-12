#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: __ ___ 2023
################################################################################
# Note

# This script extracts from dbSNP genomic reference information
# (Variant Call Format; VCF) the reference SNP cluster identifier (rsID) of each
# SNP and introduces these to a set of genomic features such as summary
# statistics from genome-wide association study (GWAS). The genomic coordinates
# for the SNP in dbSNP reference and the genomic coordinates in the set of
# genomic features such as GWAS summary statistics need to use the same
# assembly of the human genome (GRCh37 or GRCh38).

# TODO: TCW; 12 March 2023
# TODO: This script is not yet complete, but it will follow the general pattern
# of script "impute_gwas_allele_frequency.sh".

################################################################################




#
