#!/bin/bash


################################################################################
# Note

# Access a single reference file in FASTA format for the Genome Reference
# Consortium (GRC) entire human genome assembly GRCh37.

# Release GRCh37.p13 was the latest release of human genome assembly GRCh37.
# The release date of GRCh37.p13 was 28 June 2013.

# Gencode is the primary source for assemblies of the human genome.
# Gencode: "https://www.gencodegenes.org/human/releases.html"

# The human genome assembly GRCh37 is also available from other sources.
# National Center for Biotechnology Information (NCBI)
# NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/
# NCBI: https://www.ncbi.nlm.nih.gov/search/all/?term=grch37.p13
# Ensembl: https://uswest.ensembl.org/index.html

################################################################################

################################################################################
# Organize arguments.
human_genome_assembly=${1} # human genome assembly, either 'grch37' or 'grch38'
path_genome_parent_container=${2} # full path to parent directory for human genome sequence

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

################################################################################
# Organize directories.
# Access information from Gencode.

if [[ "$human_genome_assembly" == "grch37" ]]; then
  echo "Accessing sequence files for human genome assembly GRCh37."
  path_genome_container="${path_genome_parent_container}/grch37"
  rm -r $path_genome_container
  mkdir -p "${path_genome_container}"
  cd $path_genome_container
  # Genome Reference Consortium (GRC) human assembly: GRCh37.p13
  # RefSeq accession: GCF_000001405.25
  # Release date: December 2013
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz # 767 Megabytes
elif [[ "$human_genome_assembly" == "grch38" ]]; then
  echo "Accessing sequence files for human genome assembly GRCh38."
  path_genome_container="${path_genome_parent_container}/grch38"
  rm -r $path_genome_container
  mkdir -p "${path_genome_container}"
  cd $path_genome_container
  # Genome Reference Consortium (GRC) human assembly: GRCh38.p13
  # RefSeq accession:
  # Release: 42
  # Release date: October 2022


  # Genome Reference Consortium (GRC) human assembly: GRCh38.p13
  # RefSeq accession: GCF_000001405.39
  # Release date: December 2021
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz # 848 Megabytes
else
  echo "invalid specification of human genome assembly"
fi



#
