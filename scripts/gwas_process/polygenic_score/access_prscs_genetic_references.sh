#!/bin/bash


################################################################################
# Notes:
# Access genetic reference files for PRS-CS and PRS-CSX.
# https://github.com/getian107/PRScs
# https://github.com/getian107/PRScsx

# PRS-CS and PRS-CSX need the following genetic references.
# 1. Linkage Disequilibrium (LD) reference panels
# 2. Single Nucleotide Polymorphism (SNP) information
# For use in PRS-CS and PRS-CSX, the creators derived these genetic references
# from the UK Biobank and from the 1000 Genomes Project. It is convenient to
# use either set of references.

# Within the parent directory container, create two child directories.
# Within first child directory, access and organize files for genetic references
# from 1000 Genomes Project.
# Within second child directory, access and organize files for genetic
# references from UK Biobank.

################################################################################

################################################################################
# Organize arguments.
path_genetic_reference=${1} # full path to directory for genetic references

###########################################################################
# Organize paths.

path_1000_genomes="$path_genetic_reference/1000_genomes"
path_uk_biobank="$path_genetic_reference/uk_biobank"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

###########################################################################
# Organize directories.

rm -r $path_genetic_reference

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_genetic_reference ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_genetic_reference
    mkdir -p $path_1000_genomes
    mkdir -p $path_uk_biobank
fi

###########################################################################
# Access references for PRS-CS and PRS-CSX.

##########
# 1000 Genomes Project
cd $path_1000_genomes

# Linkage disequilibrium (LD) reference panels.
wget https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz
wget https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz
wget https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz
wget https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz
wget https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz
tar -zxvf ldblk_1kg_afr.tar.gz
tar -zxvf ldblk_1kg_amr.tar.gz
tar -zxvf ldblk_1kg_eas.tar.gz
tar -zxvf ldblk_1kg_eur.tar.gz
tar -zxvf ldblk_1kg_sas.tar.gz
# Single Nucleotide Polymorphism (SNP) information.
wget https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3

##########
# UK Biobank
cd $path_uk_biobank

# Linkage disequilibrium (LD) reference panels.
wget https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz
wget https://www.dropbox.com/s/y7ruj364buprkl6/ldblk_ukbb_amr.tar.gz
wget https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz
wget https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz
wget https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz
tar -zxvf ldblk_ukbb_afr.tar.gz
tar -zxvf ldblk_ukbb_amr.tar.gz
tar -zxvf ldblk_ukbb_eas.tar.gz
tar -zxvf ldblk_ukbb_eur.tar.gz
tar -zxvf ldblk_ukbb_sas.tar.gz
# Single Nucleotide Polymorphism (SNP) information.
wget https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3



#
