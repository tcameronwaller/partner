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
path_genotype_source_vcf=${1} # full path to source genotype file in VCF format
path_genotype_product_bim_container=${2} # full path to directory for product genotype files in BIM format
name_file_product_bim=${3} # name for product file in BIM format
report=${4} # whether to print reports

###########################################################################
# Organize paths.

# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_plink2=$(<"./tools_plink2.txt")

################################################################################
# Read information from VCF format file to PLINK2 and write file in BIM format.
# https://www.cog-genomics.org/plink/2.0/data#make_pgen
# Include the 'zs' modifier for '--make-just-bim' to apply Z-standard compression.

cd $path_genotype_product_bim_container

$path_plink2 \
--memory 90000 \
--threads 4 \
--vcf $path_genotype_source_vcf \
--make-just-bim \
--out "${name_file_product_bim}"

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "convert_vcf_to_plink_bim.sh"
  echo "----------"
  echo "___ report:"
  #cat $path_heritability_report_suffix
  echo "----------"
fi
