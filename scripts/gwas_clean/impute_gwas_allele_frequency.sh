#!/bin/bash

################################################################################
# Notes:

##########
# Tool: GWAS2VCF
# PubMed: 33441155
# GWAS-VCF format specification: https://github.com/MRCIEU/gwas-vcf-specification
# Host of GWAS2VCF: https://github.com/MRCIEU/gwas2vcf
# Documentation for GWAS2VCF: https://mrcieu.github.io/gwas2vcf
# Host of GWASGlue: https://github.com/MRCIEU/gwasglue
# Examples of analyses: https://mrcieu.github.io/gwasglue/articles/

# Last execution: TCW; 17 February 2023

################################################################################

################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

path_directory_sequence="${path_directory_parent}/genome_sequence"
#path_directory_dbsnp="${path_directory_parent}/dbsnp"

# Initialize directory.
#rm -r $path_directory_sequence
#rm -r $path_directory_dbsnp
mkdir -p $path_directory_sequence
#mkdir -p $path_directory_dbsnp

###########################################################################
# Execute procedure.

# View header of 1000 Genomes file in Variant Call Format (VCF) to determine
# available tags for extraction.
$path_bcftools head $path_file_1000_genomes_vcf

# View records within 1000 Genomes file.
$path_bcftools view --no-header --regions 7:1000-1010 $path_file_1000_genomes_vcf

# Filter to sites (loci) that only have two allelic variants (biallelic sites).
#$path_bcftools norm --multiallelics +snps $path_file_1000_genomes_vcf | $path_bcftools view --no-header --min-alleles 2 --max-alleles 2 --types snps | head -10
$path_bcftools norm --multiallelics +snps $path_file_1000_genomes_vcf | $path_bcftools view --min-alleles 2 --max-alleles 2 --types snps > $path_ftemp_biallelic_sites

# Extract relevant information to a flat text table.
$path_bcftools query -f '%ID %CHROM %POS %REF %ALT %INFO/EUR\n' $path_ftemp_biallelic_loci | head -100

# TODO: TCW; 21 February 2023
# TODO: organize steps above in an executable script that is conveniently callable for a single set of GWAS sum stats
# TODO: - - Include turn-on-able reports that show the first few lines of the flat text table from the bcftools query, for example.
# TODO: use awk to construct a unique row identifier with rsID_chrom_position_alternate-allele
# TODO: use that unique row identifier to merge in awk to the GWAS summ stats
# TODO: adjust format appropriately


#
