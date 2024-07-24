#!/bin/bash

#chmod u+x script.sh

# Echo each command to console.
set -x

# Set working directory.
cd ~



##################################################
# Install tools for biomedical data science
##################################################



##################################################
# General

##########
# name: Git, GitHub
# site: https://github.com/
# date, installation:
# date, review: 19 June 2024
# After renaming a repository, redirect git commits to new repository on GitHub.
git remote set-url origin https://github.com/tcameronwaller/partner.git
git remote set-url origin https://github.com/tcameronwaller/psychiatry_biomarkers.git

##########
# name: R Project for Statistical Computing
# site: https://cran.r-project.org/
# date, installation:
# date, review:
# Install R from the CRAN Repository (not the Ubuntu repository).
sudo apt update -qq
sudo apt install --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt install --no-install-recommends r-base
sudo add-apt-repository ppa:c2d4u.team/c2d4u4.0+

##########
# name: Python
# site: https://www.python.org/
# date, installation:
# date, review:
# note: a separate file describes installation of Python and its packages

##########
# name: Cytoscape
# site: https://cytoscape.org/index.html
# date, installation: 22 February 2024
# date, review: 22 February 2024
# 1. Install Java 17
apt update
apt upgrade
apt install openjdk-17-jdk openjdk-17-jre
java -version
# 2. Install Cytoscape
cd ~/Downloads
wget https://github.com/cytoscape/cytoscape/releases/download/3.10.1/Cytoscape_3_10_1_unix.sh
bash ~/Downloads/Cytoscape_3_10_1_unix.sh
# installation path: /home/tcameronwaller/Cytoscape_v3.10.1



##################################################
# Transcript

##########
# name: STAR
# site: ___
# documentation: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# GitHub: https://github.com/alexdobin/STAR
# version, latest: v2.7.11b
# - date, release: 25 January 2024
# installation:
# - date, installation: 20 June 2024
# - version, installation: v2.7.11b
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz
tar -xzvf ./2.7.11b.tar.gz -C ./
rm ./2.7.11b.tar.gz
cd ./STAR-2.7.11b/source
make
make STAR # "make: 'STAR' is up to date."
# Execution from compilation of source.
/.../tool/STAR-2.7.11b/source/STAR --help
# Execution from pre-compiled version.
/.../tool/STAR-2.7.11b/bin/Linux_x86_64/STAR --help
# Execution from pre-compiled version with static executables without external dependencies.
/.../tool/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --help


##################################################
# Genotype


##########
# name: SamTools
# site: http://www.htslib.org/
# documentation: http://www.htslib.org/doc/#manual-pages
# GitHub: https://github.com/samtools/samtools
# version, latest: v1.20
# - date, release: 15 April 2024
# installation:
# - system: NCSA, mForge, endocrinology workspace
# - date, installation: 1 July 2024
# - version, installation: v1.20
cd /.../tool
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -xjvf ./samtools-1.20.tar.bz2 -C ./
rm ./samtools-1.20.tar.bz2
cd ./samtools-1.20
pwd
./configure --prefix="/.../tool/samtools-1.20" # requires full, absolute path to directory
make
make install
/.../tool/samtools-1.20/bin/samtools --help # must execute with path


##########
# name: BCFTools
# site: http://www.htslib.org/
# documentation: http://www.htslib.org/doc/#manual-pages
# GitHub: http://github.com/samtools/bcftools
# version, latest: v1.20
# - date, release: 15 April 2024
# installation:
# - system: NCSA, mForge, endocrinology workspace
# - date, installation: 1 July 2024
# - version, installation: v1.20
cd /.../tool
wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
tar -xjvf ./bcftools-1.20.tar.bz2 -C ./
rm ./bcftools-1.20.tar.bz2
cd ./bcftools-1.20
pwd
./configure --prefix="/.../tool/bcftools-1.20" # requires full, absolute path to directory
make
make install
/.../tool/bcftools-1.20/bin/bcftools --help # must execute with path


##########
# name: HTSLib (BGZip, Tabix, etc)
# site: http://www.htslib.org/
# documentation: http://www.htslib.org/doc/#manual-pages
# GitHub: https://github.com/samtools/htslib
# version, latest: v1.20
# - date, release: 15 April 2024
# installation:
# - system: NCSA, mForge, endocrinology workspace
# - date, installation: 1 July 2024
# - version, installation: v1.20
cd /.../tool
wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
tar -xjvf ./htslib-1.20.tar.bz2 -C ./
rm ./htslib-1.20.tar.bz2
cd ./htslib-1.20
pwd
./configure --prefix="/.../tool/htslib-1.20" # requires full, absolute path to directory
make
make install
/.../tool/htslib-1.20/bin/bgzip --help # must execute with path
/.../tool/htslib-1.20/bin/tabix --help # must execute with path


##########
# name: GCTB
# site: https://cnsgenomics.com/software/gctb/download/README.html
# date, installation: 13 January 2023
# date, review:
# A package of tools for Genome-wide Complex Trait Bayesian (GCTB) analysis
# Includes the tool SBayesR for calculation of Polygenic Scores (PGS).
# Authors provided a "statically linked 64-bit Linux executable, gctb".
# Installation instructions for custom, local compilation (https://cnsgenomics.com/software/gctb/download/README.html).
cd ./ # navigate to directory in which to install program.
wget "https://cnsgenomics.com/software/gctb/download/gctb_2.04.3_Linux.zip"
unzip "./gctb_2.04.3_Linux.zip"
./gctb_2.04.3_Linux/gctb --help


##########
# name: LDSC
# site: https://github.com/bulik/ldsc
# date, installation:
# date, review:
# A package of tools to estimate autosome-wide (excluding sex chromosomes) SNP
# heritability and genetic correlation by Linkage Disequilibrium (LD) Score
# Regression.
# PubMed: ___
# documentation:
# LDSC is a Python package, and it is necessary to run this program within a
# special Python environment.
# Refer to the script "install_python_virtual_environments_packages.sh".
# Navigate to the directory in which to install program.
cd ./ldsc/
# Copy GitHub repository.
git clone https://github.com/bulik/ldsc.git

##########
# name: PRS-CSX
# site: https://github.com/getian107/PRScsx
# date, installation:
# date, review:
# Install PRS-CS and PRS-CSX.
cd ./prs_cs/ # Navigate to the directory in which to install program.
git clone https://github.com/getian107/PRScs.git
git clone https://github.com/getian107/PRScsx.git

##########
# name: GWAS2VCF
# site: https://github.com/MRCIEU/gwas2vcf
# date, installation: 23 January 2023
# date, review:
# A package of tools to translate GWAS summary statistics to the GWAS-VCF
# format.
# PubMed: 33441155
# GWAS-VCF format specification: https://github.com/MRCIEU/gwas-vcf-specification
# Host of GWAS2VCF: https://github.com/MRCIEU/gwas2vcf
# Documentation for GWAS2VCF: https://mrcieu.github.io/gwas2vcf/install/#dbsnp
# GWAS2VCF is a Python package, and it is necessary to run this program within a
# special Python environment.
# Refer to the script "install_python_virtual_environments_packages.sh".
# Access the GWAS2VCF package.
cd ./tools
mkdir -p ./gwas2vcf
cd ./gwas2vcf # Navigate to the directory in which to install program.
git clone https://github.com/MRCIEU/gwas2vcf.git
python3 ./gwas2vcf/main.py # Execute within a Python virtual environment with dependencies.

##########
# name: GZ-Sort
# site: https://github.com/keenerd/gz-sort
# date, installation: 7 February 2023
# date, review:
# site: http://kmkeen.com/gz-sort/
# repository: https://github.com/keenerd/gz-sort
git clone https://github.com/keenerd/gz-sort; cd gz-sort; make; ./gz-sort -h

##########
# name: CrossMap
# site: https://github.com/liguowang/CrossMap
# date, installation:
# date, review:
# CrossMap is a Python package, and it is necessary to run this program within a
# special Python environment.
# Refer to the script "install_python_virtual_environments_packages.sh".



#
