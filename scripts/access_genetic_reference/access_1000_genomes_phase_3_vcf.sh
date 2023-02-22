#!/bin/bash


################################################################################
# Notes:

# This script accesses genomic variation information in Variant Call Format
# (vcf) from Phase 3 of the 1000 Genomes project.
# EnsEMBL hosts these data.
# This variant information uses the GRCh37 assembly of the human genome.

# EnsEMBL host website.
# Site: https://grch37.ensembl.org/info/data/ftp/index.html

# Ensembl File Transfer Protocol (FTP) server.
# Site: https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/

# Review: TCW; __ February 2023
# Last Accession: __ February 2023 (TCW; to NCSA server)

################################################################################



################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

# Initialize directory.
mkdir -p $path_directory_parent

cd $path_directory_parent

###########################################################################
# Execute procedure.

# Echo each command to console.
#set -x
# Suppress echo each command to console.
#set +x

###########################################################################
# Organize directories.
# Access reference information from NCBI dbSNP.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Accessing 1000 Genomes genomic variation information from EnsEMBL."
  echo "Human Genome Assembly GRCh37."
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Access the specific files and save within the directory.
wget "https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/README" --directory-prefix $path_directory_parent --content-disposition --no-check-certificate --show-progress
wget "https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz" --directory-prefix $path_directory_parent --content-disposition --no-check-certificate --show-progress
wget "https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi" --directory-prefix $path_directory_parent --content-disposition --no-check-certificate --show-progress


################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "access_1000_genomes_phase_3_vcf.sh"
  echo "----------"
fi





#
