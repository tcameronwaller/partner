#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 10 March 2023
################################################################################
# Notes

# Source Format
# Description: Format of SNP effect weights from SBayesR.
# Documentation site: https://cnsgenomics.com/software/gctb/#Tutorial
# File suffix: ".snpRes"
# File type: text
# File compression: none
# Delimiter: Tab
# Chromosome base position coordinate system: base 1
#   Site: https://www.biostars.org/p/84686/
#   Note: Coordinates designate 1-based integer position of each base
# Columns: Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP (TCW; 2023-03-10)
#          1  2    3     4        5  6  7     8        9  10

# Product Format
# Description: Browser Extensible Data (BED) format for UCSC Genome Browser and CrossMap
# Documentation site: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
# Documentation site: https://crossmap.sourceforge.net/#convert-bed-format-files
# File suffix: ".bed.gz"
# File type: text
# File compression: Gzip
# Delimiter: Tab
# Chromosome base position coordinate system: base 0
#   Site: https://www.biostars.org/p/84686/
#   Note: Coordinates designate 0-based integer range flanking base or range of bases
# Columns: chrom chromStart chromEnd Name A1 A2 A1Frq A1Effect SE PIP
#          1     2          3        4    5  6  7     8        9  10

# review: TCW; 10 March 2023

###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize arguments.

path_file_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
report=${3} # whether to print reports

################################################################################
# Organize paths.

name_base_file_product="$(basename $path_file_product ".bed.gz")"
path_directory_product="$(dirname $path_file_product)"
path_directory_product_temporary="${path_directory_product}/temporary_translation_${name_base_file_product}" # hopefully unique

path_file_temporary_format="${path_directory_product_temporary}/${name_base_file_product}_temporary.bed"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product_temporary

# Remove any previous version of the product file.
rm $path_file_product

###########################################################################
# Execute procedure.

##########
# Translate format of information about genomic features.
# AWK interprets a single space delimiter (FS=" ") as any white space.

echo -e "chrom\tchromStart\tchromEnd\tName\tA1\tA2\tA1Frq\tA1Effect\tSE\tPIP" > $path_file_temporary_format
cat $path_file_source | awk 'BEGIN { FS=" "; OFS="\t" } NR > 1 {
  print $3, ($4 - 1), ($4), $2, toupper($5), toupper($6), $7, $8, $9, $10
}' >> $path_file_temporary_format

# Compress file format.
gzip -cvf $path_file_temporary_format > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Translate format of SNP effect weights from SBayesR to"
  echo "BED format for UCSC Genome Browser and CrossMap."
  echo "path to source SNP effect weights file: " $path_file_source
  echo "path to product SNP effect weights file: " $path_file_product
  echo "table before format translation:"
  head -10 $path_file_source
  echo "table after format translation:"
  head -10 $path_file_temporary_format
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
