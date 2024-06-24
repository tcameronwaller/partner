#!/bin/bash

################################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, first execution: 6 February 2023
# Date, last execution: 20 June 2024
# Review: TCW; 20 June 2024
################################################################################
# Note

# TODO: TCW; 24 June 2024
# This script is incomplete.
# Update the "wget" accession paths.


# Access a single file in FASTA format for the assemblies GRCh37 and GRCh38 of
# the human reference genome.

# The source for these accessions is GENCODE (https://www.gencodegenes.org/).

# Genome Reference Consortium human genome assembly 38 (GRCh38)
# - version, latest: GRCh38 patch 14 (GRCh38.p14), GENCODE release 46
# - date, release original GRCh38.p14: 3 February 2022
# - date, GENCODE release 46: 1 May 2024
# - format: FASTA

# Genome Reference Consortium human genome assembly 37 (GRCh37)
# - version, latest: GRCh37 patch 13 (GRCh38.p13), GENCODE release 19
# - date, release original GRCh37.p13: 28 June 2013
# - date, GENCODE release 19: 1 December 2013
# - format: FASTA

# Reference
# - site, Genome Reference Consortium: https://www.ncbi.nlm.nih.gov/grc
# - site, Gencode: https://www.gencodegenes.org/human/
# - site, Ensembl: https://useast.ensembl.org/Homo_sapiens/Info/Annotation

# Gencode is the primary source for assemblies of the human genome.
# Gencode: "https://www.gencodegenes.org/human/releases.html"

# The human genome assembly GRCh37 is also available from other sources.
# National Center for Biotechnology Information (NCBI)
# NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/
# NCBI: https://www.ncbi.nlm.nih.gov/search/all/?term=grch37.p13
# Ensembl: https://uswest.ensembl.org/index.html

# Review: TCW; 20 June 2024

################################################################################



################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

path_directory_grch37="${path_directory_parent}/gencode/grch37"
path_directory_grch38="${path_directory_parent}/gencode/grch38"

# Initialize directory.
rm -r $path_directory_grch37
rm -r $path_directory_grch38
mkdir -p $path_directory_grch37
mkdir -p $path_directory_grch38

###########################################################################
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
  echo "Accessing sequence reference files for human genome assembly GRCh38."
  echo "----------"
  echo "----------"
  echo "----------"
fi
# Accession.
cd $path_directory_grch38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.p13.genome.fa.gz # 849 Megabytes






##########
# GRCh37

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Accessing sequence reference files for human genome assembly GRCh37."
  echo "----------"
  echo "----------"
  echo "----------"
fi
cd $path_directory_grch37
# Genome Reference Consortium (GRC) human assembly: GRCh37.p13
# Description: https://www.gencodegenes.org/human/release_19.html
# RefSeq accession: GCF_000001405.25
# Release date: December 2013
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz # 768 Megabytes



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "access_human_genome_sequence_gencode_grch37_grch38_fasta.sh"
  echo "----------"
fi



#
