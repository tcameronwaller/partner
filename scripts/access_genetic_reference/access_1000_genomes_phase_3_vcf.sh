#!/bin/bash

################################################################################
# Notes:

# This script accesses genomic variation information in Variant Call Format
# (VCF) from an EnsEMBL host version of Phase 3 of the 1000 Genomes project.
# This genomic variation information uses the GRCh37 assembly of the human
# genome.
# This genomic variation information includes frequencies of alternate alleles
# calculated in superpopulations of the 1000 Genomes Project, such as the
# European superpopulation.

# EnsEMBL host website.
# Site: https://grch37.ensembl.org/info/data/ftp/index.html

# EnsEMBL File Transfer Protocol (FTP) server.
# Site: https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/

# Review: TCW; 22 February 2023
# Last Accession: TCW; 22 February 2023 (to NCSA mForge server)

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
  echo "Accessing EnsEMBL host version of genomic variation information in"
  echo "Phase 3 of the 1000 Genomes Project."
  echo "Human Genome Assembly: GRCh37"
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Access the specific files and save within the directory.

# File date: 2021-03-20; File size: 5.5 K
wget --directory-prefix $path_directory_parent --content-disposition --no-check-certificate "https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/README"

# File date: 2021-04-04; File size: 1.6 G
wget --directory-prefix $path_directory_parent --content-disposition --no-check-certificate "https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"

# File date: 2021-04-04; File size: 1.7 M
wget --directory-prefix $path_directory_parent --content-disposition --no-check-certificate "https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi"


################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "access_1000_genomes_phase_3_vcf.sh"
  echo "----------"
fi





#
