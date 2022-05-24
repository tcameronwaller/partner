#!/bin/bash


################################################################################
# Note

# Access a single reference file in FASTA format for the Genome Reference
# Consortium (GRC) entire human genome assembly GRCh37.

# Release GRCh37.p13 was the latest release of human genome assembly GRCh37.
# The release date of GRCh37.p13 was 28 June 2013.

# Gencode is the primary source for assemblies of the human genome.
# Gencode: "https://www.gencodegenes.org/human/release_19.html"

# The human genome assembly GRCh37 is also available from other sources.
# National Center for Biotechnology Information (NCBI)
# NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/
# NCBI: https://www.ncbi.nlm.nih.gov/search/all/?term=grch37.p13
# Ensembl: https://uswest.ensembl.org/index.html

################################################################################

# TODO: TCW; 23 May 2022
# TODO: the access script feels too unruly to me with both GRCh37 and GRCh38...
# TODO: maybe pass a flag for the option to switch between the two?


################################################################################
# Organize arguments.
path_human_genome_grch37=${1} # full path to directory for human genome assembly GRCh37

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

###########################################################################
# Organize directories.

rm -r $path_human_genome_grch37

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_human_genome_grch37 ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_human_genome_grch37
fi

cd $path_human_genome_grch37

###########################################################################
# Access reference information from NCBI dbSNP.

# Human genome assembly: GRCh37.p13
# Release date: December 2013

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz # 767 Megabytes



#
