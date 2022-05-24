#!/bin/bash


################################################################################
# Note

# Chain files specify the translations of coordinates for chromosomes and base
# pair positions of Single Nucleotide Polymorphisms (SNPs) from a source human
# genome assembly to a target human genome assembly.

# The CrossMap tool for translations between human genome assemblies can use
# chain files from the UCSC Genome Browser or Ensembl (among others).
# Documentation for Python 2 CrossMap: "https://pythonhosted.org/CrossMap/"
# Documentation for Python 3 CrossMap: "https://sourceforge.net/projects/crossmap/"
# Documentation for Python 3 CrossMap: "http://crossmap.sourceforge.net/"

# Genome Reference Consortium (GRC) human assembly: GRCh37.p13
# Genome Reference Consortium (GRC) human assembly: GRCh38.p14

# File Transfer Protocol (FTP)
# http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/

################################################################################

# TODO: TCW; 24 May 2022
# TODO: also download the "hg19" to/from "hg38" chain files from UCSC
# TODO: links in the CrossMap documentation

################################################################################
# Organize arguments.
path_assembly_chain_container=${1} # full path to directory for chain files

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

###########################################################################
# Organize directories.
# Access chain files to map between assemblies of the human genome.

rm -r $path_assembly_chain_container
mkdir -p "${path_assembly_chain_container}"
cd $path_assembly_chain_container

# Source human genome assembly: GRCh37
# Target human genome assembly: GRCh38
# File date: 25 July 2014
# Host: Ensembl
wget http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz # 279 Kilobytes

# Source human genome assembly: GRCh38
# Target human genome assembly: GRCh37
# File date: 25 July 2014
# Host: Ensembl
wget http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz # 713 Kilobytes




#
