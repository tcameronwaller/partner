#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 6 February 2023
# Date, last execution or modification: 12 July 2024
# Review: TCW; 12 July 2024
###############################################################################
# Note

###############################################################################
# Organize arguments.

path_directory=${1} # full path to directory within which to save files
report=${2} # whether to print reports

###############################################################################
# Organize paths.

# Initialize directory.
mkdir -p $path_directory

###############################################################################
# Execute procedure.



##########
# GRCh38
# Genome Reference Consortium human genome assembly 38 (GRCh38)
# - version, latest: GRCh38 patch 14 (GRCh38.p14), GENCODE release 46
# - date, release original GRCh38.p14: 3 February 2022
# - date, GENCODE release 46: 1 May 2024
# - format: FASTA
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Accessing annotation reference files for human genome assembly GRCh38."
  echo "----------"
  echo "----------"
  echo "----------"
fi
# Accession.
cd $path_directory
# Primary assembly.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz # 49 Megabytes
# Complete assembly.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz # 53 Megabytes



###############################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "access_human_genome_annotation_gencode_grch38_gtf.sh"
  echo "----------"
fi



###############################################################################
# End.
