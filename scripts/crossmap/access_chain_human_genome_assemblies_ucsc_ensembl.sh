#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 24 May 2022
# Date, last execution: 24 May 2023
# Review: TCW; 24 May 2023
################################################################################
# Note

# This script accesses chain files for the translation of genomic coordinates
# from one assembly of the human genome to another.
# Chain files specify the translations of coordinates for chromosomes and base
# pair positions of Single Nucleotide Polymorphisms (SNPs) from a source human
# genome assembly to a target human genome assembly.

# University of California Santa Cruz (UCSC): https://genome.ucsc.edu/

# The CrossMap tool for translations between human genome assemblies can use
# chain files from the UCSC Genome Browser or Ensembl (among others).
# Documentation for Python 2 CrossMap: "https://pythonhosted.org/CrossMap/"
# Documentation for Python 3 CrossMap: "https://sourceforge.net/projects/crossmap/"
# Documentation for Python 3 CrossMap: "http://crossmap.sourceforge.net/"

# Genome Reference Consortium (GRC) human assembly: GRCh37.p13
# Genome Reference Consortium (GRC) human assembly: GRCh38.p14

################################################################################



################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

path_directory_ucsc="${path_directory_parent}/ucsc"
path_directory_ensembl="${path_directory_parent}/ensembl"

# Initialize directory.
rm -r $path_directory_ucsc
rm -r $path_directory_ensembl
mkdir -p $path_directory_ucsc
mkdir -p $path_directory_ensembl

###########################################################################
# Execute procedure.

# Echo each command to console.
#set -x
# Suppress echo each command to console.
#set +x



##########
# Assembly chain files from UCSC.
# http://hgdownload.soe.ucsc.edu/goldenPath/
# Source human genome assembly: NCBI36 (hg18)
# Target human genome assembly: GRCh37 (hg19)
# File date: 26 July 2010
# File size: 137 kilobytes
# Host: UCSC
wget --directory-prefix $path_directory_ucsc --content-disposition --no-check-certificate "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz"
# Source human genome assembly: NCBI36 (hg18)
# Target human genome assembly: GRCh38 (hg38)
# File date: 19 February 2014
# File size: 336 kilobytes
# Host: UCSC
wget --directory-prefix $path_directory_ucsc --content-disposition --no-check-certificate "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz"
# Source human genome assembly: GRCh37 (hg19)
# Target human genome assembly: GRCh38 (hg38)
# File date: 31 December 2013
# File size: 222 kilobytes
# Host: UCSC
wget --directory-prefix $path_directory_ucsc --content-disposition --no-check-certificate "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
# Source human genome assembly: GRCh38 (hg38)
# Target human genome assembly: GRCh37 (hg19)
# File date: 31 December 2013
# File size: 1.2 Megabytes
# Host: UCSC
wget --directory-prefix $path_directory_ucsc --content-disposition --no-check-certificate "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"



##########
# Assembly chain files from Ensembl.
# http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/
# Source human genome assembly: NCBI36 (hg18)
# Target human genome assembly: GRCh37 (hg19)
# File date: 28 July 2014
# File size: 30 kilobytes
# Host: Ensembl
wget --directory-prefix $path_directory_ensembl --content-disposition --no-check-certificate "http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/NCBI36_to_GRCh37.chain.gz"
# Source human genome assembly: NCBI36 (hg18)
# Target human genome assembly: GRCh38 (hg38)
# File date: 25 July 2014
# File size: 174 kilobytes
# Host: Ensembl
wget --directory-prefix $path_directory_ensembl --content-disposition --no-check-certificate "http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/NCBI36_to_GRCh38.chain.gz"
# Source human genome assembly: GRCh37
# Target human genome assembly: GRCh38
# File date: 25 July 2014
# File size: 279 kilobytes
# Host: Ensembl
wget --directory-prefix $path_directory_ensembl --content-disposition --no-check-certificate "http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz"
# Source human genome assembly: GRCh38
# Target human genome assembly: GRCh37
# File date: 25 July 2014
# File size: 713 kilobytes
# Host: Ensembl
wget --directory-prefix $path_directory_ensembl --content-disposition --no-check-certificate "http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz"



#
