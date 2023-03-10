#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 10 March 2023
################################################################################
# Note

# Documentation for Python 2 CrossMap: "https://pythonhosted.org/CrossMap/"
# Documentation for Python 3 CrossMap: "https://sourceforge.net/projects/crossmap/"
# Documentation for Python 3 CrossMap: "http://crossmap.sourceforge.net/"




# The CrossMap tool translates coordinates for chromosomes, base-pair positions,
# and alleles of Single Nucleotide Polymorphisms (SNPs) from a source human
# genome assembly to a target human genome assembly.
# The Crossmap tool does not change annotations for SNPs (such as identifiers
# from dbSNP).

# For translations of genotype files in Variant Call Format (VCF), CrossMap
# requires a single reference file in FASTA format for the entire target human
# genome assembly.

# Information on how to access a single FASTA file for assembly GRCh37 of the
# Human Genome: "https://www.biostars.org/p/338914/#339258"

# In the process of translating coordinates between genome assemblies, some
# genetic features (SNPs etc) change chromosomes. For this reason it might be
# appropriate to combine VCF files for all chromosomes (using BCFTools "concat"
# command) before the translation and then split by chromosome after the
# translation (using BCFTools "view" command with "--regions" option).

################################################################################

################################################################################
# Organize arguments.
path_file_source=${1} # full path to source file in UCSC Browser Extensible Data (BED) format
path_file_product=${2} # full path to product file in UCSC Browser Extensible Data (BED) format
path_file_assembly_translation_chain=${3} # full path to chain file for assembly translation
path_file_genome_sequence=${4} # full path to file in FASTA format without compression for genome sequence that matches product assembly
threads=${5} # count of processing threads to use
report=${8} # whether to print reports

################################################################################
# Organize paths.

cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_environment_crossmap="${path_waller_tools}/python/environments/crossmap"

name_base_file_product="$(basename $path_file_product .bed.gz)"
path_directory_product="$(dirname $path_file_product)"

# Initialize directory.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Execute procedure.

# Regulate concurrent or parallel process threads on node cores.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads


# Translate coordinates for chromosomes and base pair positions between human
# genome assemblies.

# Activate Virtual Environment.
source "${path_environment_crossmap}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s

# Read VCF file without compression in CrossMap.
# Translate coordinates between genome assemblies.
# Write VCF file without compression.
# CrossMap uses GZip compression ("--compress" command).
# Do not use GZip compression in order to simplify reading into BCFTools.
CrossMap.py \
bed \
--chromid a \
$path_file_assembly_translation_chain \
$path_file_source \
$path_file_product

# Deactivate Virtual Environment.
deactivate
which python3


# Remove temporary, intermediate files.
#rm -r $path_directory_product_temporary



#
