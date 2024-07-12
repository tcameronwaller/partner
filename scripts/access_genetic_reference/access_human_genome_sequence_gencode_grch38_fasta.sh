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

# Access a single file in FASTA format for the reference sequence of assembly
# GRCh38 of the human genome from GENCODE.

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
# - site, GenBank: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
# - site, Gencode: https://www.gencodegenes.org/human/
# - site, Ensembl: https://useast.ensembl.org/Homo_sapiens/Info/Annotation

# Gencode is a source for sequence and annotation of assemblies of the human
# genome.
# Gencode: "https://www.gencodegenes.org/human/releases.html"

# The human genome assembly GRCh37 is also available from other sources.
# National Center for Biotechnology Information (NCBI)
# NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/
# NCBI: https://www.ncbi.nlm.nih.gov/search/all/?term=grch37.p13
# Ensembl: https://uswest.ensembl.org/index.html

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
  echo "Accessing sequence reference files for human genome assembly GRCh38."
  echo "----------"
  echo "----------"
  echo "----------"
fi
# Accession.
cd $path_directory
# Primary assembly.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz # 807 Megabytes
# Complete assembly.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz # 856 Megabytes



###############################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo $0 # Print full file path to script.
  echo "access_human_genome_sequence_gencode_grch38_fasta.sh"
  echo "----------"
fi



###############################################################################
# End.
