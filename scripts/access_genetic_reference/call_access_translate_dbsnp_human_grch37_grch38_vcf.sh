#!/bin/bash


################################################################################
# Notes:

# Last execution: __ February 2023 (TCW; on NCSA server)

# Review: TCW; 6 February 2023

################################################################################

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_parent="${path_directory_reference}/dbsnp"
path_directory_process=$(<"./process_psychiatric_metabolism.txt")

# Files.
path_file_translation_grch37="${path_directory_process}/scripts/access_genetic_reference/translations_chromosomes_refseq_grch37p13.txt"
path_file_translation_grch38="${path_directory_process}/scripts/access_genetic_reference/translations_chromosomes_refseq_grch38p14.txt"

# Scripts.
path_file_script_access="${path_directory_process}/promiscuity/scripts/access_genetic_reference/access_dbsnp_human_grch37_grch38_vcf.sh"
path_file_script_translate="${path_directory_process}/promiscuity/scripts/access_genetic_reference/translate_chromosomes_dbsnp_human_grch37_grch38_vcf.sh"
path_file_script_bcftools_annotate="${path_directory_process}/promiscuity/scripts/bcftools/translate_chromosome_identifiers_vcf.sh"

# Initialize directories.
#rm -r $path_directory_parent
mkdir -p $path_directory_parent
################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.



if true; then
  # Call script for accession.
  /usr/bin/bash $path_file_script_access \
  $path_directory_parent \
  $report
fi

if false; then
  # Call script for translation.
  /usr/bin/bash $path_file_script_translate \
  $path_directory_parent \
  $path_file_translation_grch37 \
  $path_file_translation_grch38 \
  $path_file_script_bcftools_annotate \
  $path_bcftools \
  $report
fi

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "call_access_translate_dbsnp_human_grch37_grch38_vcf.sh"
  echo "----------"
fi





#
