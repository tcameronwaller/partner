#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: __ March 2023
# Date, last execution: __ March 2023
# Review: TCW; 15 March 2023
################################################################################
# Note

# This script translates coordinates for chromosomes, base-pair positions, and
# alleles for a set of genomic features such as single-nucleotide polymorphisms
# between assemblies of the human genome.

# Documentation for Python 2 CrossMap: "https://pythonhosted.org/CrossMap/"
# Documentation for Python 3 CrossMap: "https://sourceforge.net/projects/crossmap/"
# Documentation for Python 3 CrossMap: "http://crossmap.sourceforge.net/"
# Documentation for Python 3 CrossMap: "https://crossmap.readthedocs.io/en/latest/"
# The CrossMap tool translates coordinates for chromosomes, base-pair positions,
# and alleles of Single Nucleotide Polymorphisms (SNPs) from a source human
# genome assembly to a target human genome assembly.

# The Crossmap tool does not change annotations for SNPs such as the reference
# SNP cluster identifier (rsID) from dbSNP; however, these rsIDs are semi-stable
# designations of a genomic locus that is not specific to assembly of the human
# genome (GRCh37 or GRCh38). Hence the rsID of a SNP should not need to change
# after CrossMap translates its coordinates from GRCh37 to GRCh38 assemblies. It
# might still be reasonable to check for errors in these rsIDs in the new
# assembly.

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

# The product file in BED format that CrossMap creates does not have a header.
# This script introduces a header.

################################################################################

################################################################################
# Organize arguments.
path_file_source=${1} # full path to source file in UCSC Browser Extensible Data (BED) format
path_file_product=${2} # full path to product file in UCSC Browser Extensible Data (BED) format
path_file_chain=${3} # full path to chain file for assembly translation
threads=${4} # count of processing threads to use
report=${5} # whether to print reports

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_environment_crossmap="${path_waller_tools}/python/environments/crossmap"

path_directory_product="$(dirname $path_file_product)"
name_base_file_product="$(basename $path_file_product .bed.gz)"
path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}" # hopefully unique

# Files.
path_file_temporary_map="${path_directory_product_temporary}/${name_base_file_product}_map.bed"
path_file_temporary_map_header="${path_directory_product_temporary}/${name_base_file_product}_map_header.bed"
path_file_product_unmap="${path_directory_product}/${name_base_file_product}_unmap.txt"

# Initialize directory.
#rm -r $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_product_temporary
cd $path_directory_product_temporary

################################################################################
# Execute procedure.

# Activate Virtual Environment.
source "${path_environment_crossmap}/bin/activate"
echo "confirm Python Virtual Environment path..."
which python3
sleep 5s
# Regulate concurrent or parallel process threads on node cores.
# Force Python program (especially SciPy) not to use all available cores on a
# cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads
# Call CrossMap.
# Translate coordinates for chromosomes, base pair positions, and alleles
# between human genome assemblies.
# CrossMap uses GZip compression ("--compress" command).
# I think that CrossMap by default does not compress product BED files.
CrossMap.py \
bed \
--chromid a \
--unmap-file $path_file_product_unmap \
$path_file_chain \
$path_file_source \
$path_file_temporary_map
# Deactivate Virtual Environment.
deactivate
which python3

# Introduce the same header from the source file to the product file.
# CrossMap does not transfer the original header.
zcat $path_file_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary_map_header
cat $path_file_temporary_map >> $path_file_temporary_map_header

# Compress file format.
gzip -cvf $path_file_temporary_map_header > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate genomic features from human genome assembly"
  echo "GRCh37 to GRCh38 in CrossMap."
  echo "path to source file: " $path_file_source
  echo "path to chain file: " $path_file_chain
  echo "path to product file: " $path_file_product
  echo "----------"
  echo "table before translation:"
  echo "Count lines: "
  zcat $path_file_source | wc -l
  zcat $path_file_source | head -10
  echo "----------"
  echo "table after translation:"
  echo "Count lines: "
  zcat $path_file_product | wc -l
  zcat $path_file_product | head -10
  echo "----------"
  echo "----------"
  echo "----------"
fi



# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
