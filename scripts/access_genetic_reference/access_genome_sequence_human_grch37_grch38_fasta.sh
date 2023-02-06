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

# Review: TCW; 6 February 2023

################################################################################



################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

path_directory_grch37="${path_directory_parent}/grch37"
path_directory_grch38="${path_directory_parent}/grch38"

# Initialize directory.
rm -r $path_directory_grch37
rm -r $path_directory_grch38
mkdir -p $path_directory_grch37
mkdir -p $path_directory_grch38

###########################################################################
# Execute procedure.



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
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz # 767 Megabytes

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
cd $path_directory_grch38
# Genome Reference Consortium (GRC) human assembly: GRCh38.p13
# Description: https://www.gencodegenes.org/human/release_42.html
# RefSeq accession: GCF_000001405.42
# Release: 42
# Release date: April 2022
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.p13.genome.fa.gz # 848 Megabytes



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "access_genome_sequence_human_grch37_grch38_fasta.sh"
  echo "----------"
fi



#
