#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 March 2023
# Date, last execution: 20 March 2023
# Review: TCW; ___
################################################################################
# Note

# This script calls a method within PLINK2 to calculate a polygenic score as
# linear combination of allelic effects across SNPs within a genotype.
# Documentation: https://www.cog-genomics.org/plink/1.9/score
# Documentation: https://www.cog-genomics.org/plink/2.0/score

# The variant identifiers in the target genotypes need to match the format of
# the variant identifiers in the allelic effects. It is likely that the
# identifiers of variants in the target genotypes are not the same as the
# reference SNP cluster identifier (rsID) that appears often in GWAS summary
# statistics.

# By default, PLINK2 "--score" method attempts to handle missing genotype SNPs
# by giving them an imputed value proportional to their cohort-specific allele
# frequency. I think these allele frequencies need to be specific to the cohort
# of target genotypes. Use the PLINK2 "--read-freq" method to read a file in
# ".afreq" format produced by PLINK2 "--freq" method on the target genotypes.
# Documentation: https://www.cog-genomics.org/plink/2.0/formats#afreq

# For calculation of cumulative allelic effects across all available SNPs in the
# genome (actually autosome, chromosomes 1-22) it is important to use the sum of
# allelic effects within each chromosome and not the mean (sum of allelic
# effects divided by the count of non-missing alleles used to calculate the
# sum).

# SNP effects
# Description: Format of SNP effects for calculation of polygenic scores
# Description: This format attempts to be compatible with team standard for GWAS summary statistics.
# Documentation site: https://www.cog-genomics.org/plink/1.9/score
# Documentation site: https://www.cog-genomics.org/plink/2.0/score
# File suffix: ".txt.gz"
# File type: text
# File compression: GZip
# Delimiter: Tab
# Header: Yes
# Chromosome base position coordinate system: base 1
#   Site: https://www.biostars.org/p/84686/
#   Note: Coordinates designate 1-based integer position of each base
# Columns: SNP CHR BP A1 A2 A1AF A1EFFECT SE P
#          1   2   3  4  5  6    7        8  9

# Polygenic scores
# Description: output file of PLINK2 "score" command is in ".sscore" format
# Documentation site: https://www.cog-genomics.org/plink/2.0/formats#sscore
# File suffix: ".sscore"

# Additional references
# https://2cjenn.github.io/PRS_Pipeline/
# https://www.biostars.org/p/9508742/

##########

# Note:
# The average would effectively just be a scaled sum score.

# Important reference!!!
# https://www.biostars.org/p/362960/
# In this reference an author of PLINK2 reports how to include additional columns
# in the PLINK2 --score report, including sums ("cols=+scoresums").
# Plan to include several of the optional columns in order to have access to this
# information and consider how to handle it in combining chromosomes.

# TODO: TCW; 16 March 2023
# TODO: I don't understand why it would be a problem to combine across chromosomes
# TODO: after PLINK2 --score calculates final scores as averages of valid alleles
# TODO: used in the calculation of the score. We can't I just calculate the sum
# TODO: of this average across chromosomes?
# TODO: reference: https://www.biostars.org/p/9508742/
# TODO: PLINK2 --score does not have the "sum" or "no-sum" modifiers.

# TODO: TCW; 16 March 2023
# TODO: another option would be to try to run PLINK2 --score on genotypes across
# TODO: all chromosomes together using the PLINK files instead of VCFs.


################################################################################
################################################################################
################################################################################

# TODO: TCW; 17 August 2022
# TODO: as the PRS-CSX documentation says, "If polygenic scores are generated by
# TODO: chromosome, use the 'sum' modifier so that they can be combined into a genome-wide score."

# TODO: TCW; 18 August 2022
# TODO: This script will accept a SINGLE file of posterior effect size report from PRS-CSX.
# TODO: This script will call PLINK2 --score for this SINGLE file.
# TODO: This script will need a driver script. Follow pattern of "3_estimate_prscs_allelic_effects.sh"


################################################################################
################################################################################
################################################################################


################################################################################
# Organize arguments.
path_file_source_effects=${1} # full path to source file in standard format of allelic effects across SNPs
path_file_source_genotypes=${2} # full path to source file in Variant Call Format (VCF) of target genotypes
path_directory_product=${3} # full path to product directory for polygenic scores of allelic effects across target genotypes
name_base_file_product=${4} # base file name for product files
threads=${5} # count of processing threads to use
report=${6} # whether to print reports

###########################################################################
# Organize paths.

# Directories.
cd ~/paths
path_plink2=$(<"./tools_plink2.txt")

#path_directory_product="$(dirname $path_file_product)"
#name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}" # hopefully unique
path_file_base_product="${path_directory_product}/${name_base_file_product}" # hopefully unique

# Files.
path_file_temporary_effects="${path_directory_product_temporary}/${name_base_file_product}.txt"
#path_file_sscore="${path_directory_product}/${name_base_file_product}.sscore"

# Initialize directory.
#rm -r $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_product_temporary
cd $path_directory_product

################################################################################
# Call PLINK2 to calculate linear combination of allelic effects across SNPs in
# target genotypes.

# Call PLINK2 help menu for the score method as this menu gives more detail than
# the documentation online.
$path_plink2 --help score

# Decompress the GWAS summary statistics.
gzip -dcvf $path_file_source_effects > $path_file_temporary_effects

# PLINK2 arguments "center", "variance-standardize", and "dominant" modify the
# scale of the polygenic scores; however, wait to adjust scale until after
# combination of scores across separate chromosomes.


#$path_plink2 \
#--memory 90000 \
#--threads $threads \
#--vcf $path_file_source_genotypes \
#--xchr-model 2 \
#--score $path_file_temporary_effects 1 4 7 header no-mean-imputation ignore-dup-ids list-variants \
#--out $name_base_file_product

$path_plink2 \
--memory 90000 \
--threads $threads \
--vcf $path_file_source_genotypes \
--xchr-model 2 \
--score $path_file_temporary_effects 1 4 header no-mean-imputation ignore-dup-ids list-variants cols=+scoresums,+denom \
--score-col-nums 7 \
--out $path_file_base_product

# Compress file format.
#gzip -cvf $path_file_sscore > $path_file_product

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
