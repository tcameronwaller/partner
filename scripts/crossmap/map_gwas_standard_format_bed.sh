#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 24 May 2023
# Date, last execution: 2 August 2023
# Review: TCW; 24 May 2023
################################################################################
# Note

# This script translates the coordinates for chromosomes, base-pair positions,
# and alleles for a set of GWAS summary statistics between assemblies of the
# human genome.

# Because GWAS summary statistics in team standard format have more columns than
# CrossMap can accommodate, some special steps are necessary.

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

# The product file in BED format that CrossMap creates does not have a header.
# This script introduces a header.

# Format of source (original) and product (final) GWAS summary statistics.
# Description: Team standard format for GWAS summary statistics.
# Variables: "path_file_source", "path_file_product"
# File suffix: ".txt.gz"
# File type: text
# File compression: Gzip
# Delimiter: Space
# Chromosome base position coordinate system: base 1
#   Site: https://www.biostars.org/p/84686/
#   Note: Coordinates designate 1-based integer position of each base
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT
#          1   2   3  4  5  6    7    8  9   10  11  12   13    14

# Format of intermediate, temporary GWAS summary statistics.
# Description: Team standard format for GWAS summary statistics.
# File suffix: ".txt.gz"
# File type: text
# File compression: Gzip
# Delimiter: Space
# Chromosome base position coordinate system: base 1
#   Site: https://www.biostars.org/p/84686/
#   Note: Coordinates designate 1-based integer position of each base
# Columns: IDENTIFIER_MERGE SNP CHR_OLD BP_OLD A1 A2 A1AF BETA SE P   N   Z   INFO NCASE NCONT
#          1                     2   3       4      5  6  7    8    9  10  11  12  13   14    15

# Format of intermediate, temporary genomic coordinates for translation.
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
# Columns: chrom chromStart chromEnd IDENTIFIER_MERGE
#          1     2          3        4



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
name_base_file_product="$(basename $path_file_product .txt.gz)"
path_directory_product_temporary="${path_directory_product}/temporary_${name_base_file_product}" # hopefully unique

# Files.
path_file_temporary_merge_source="${path_directory_product_temporary}/${name_base_file_product}_source.txt"
path_file_temporary_map_source="${path_directory_product_temporary}/${name_base_file_product}_source.bed"
path_file_temporary_map_product="${path_directory_product_temporary}/${name_base_file_product}_product.bed"
path_file_temporary_map_product_header="${path_directory_product_temporary}/${name_base_file_product}_product_header.bed"
path_file_temporary_merge_product="${path_directory_product_temporary}/${name_base_file_product}_product_before_merge.txt"
path_file_temporary_product_after_merge="${path_directory_product_temporary}/${name_base_file_product}_product_after_merge.txt"
path_file_temporary_product_filter="${path_directory_product_temporary}/${name_base_file_product}_product_filter.txt"
path_file_temporary_product_format="${path_directory_product_temporary}/${name_base_file_product}_product_format.txt"
path_file_product_unmap="${path_directory_product}/${name_base_file_product}_unmap.txt"

# Initialize directory.
#rm -r $path_directory_product # caution
rm -r $path_directory_product_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_product_temporary
cd $path_directory_product_temporary

################################################################################
# Organize parameters.

verbosity="false"

################################################################################
# Execute procedure.



##########
# 1. Introduce sequential identifier to rows in source table of GWAS summary
#     statistics.
echo "IDENTIFIER_MERGE SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_merge_source
zcat $path_file_source | awk 'BEGIN {FS = " "; OFS = " "; i = 1} NR > 1 {
  print i, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14; i++
}' >> $path_file_temporary_merge_source

# Report.
if [ "$verbosity" == "true" ] && [ "$report" == "true" ]; then
  echo "----------"
  echo "Progress Step 1: Introduction of sequential identifier for merge."
  cat $path_file_temporary_merge_source | head -5
  echo "----------"
fi



##########
# 2. Extract original genomic coordinates and sequential identifier from source
#     table and organize format in Browser Extensible Data (BED) format for
#     CrossMap.
echo -e "chrom\tchromStart\tchromEnd\tIDENTIFIER_MERGE" > $path_file_temporary_map_source
cat $path_file_temporary_merge_source | awk 'BEGIN {FS = " "; OFS = "\t"} NR > 1 {
  print $3, ($4 - 1), ($4), $1
}' >> $path_file_temporary_map_source

# Report.
if [ "$verbosity" == "true" ] && [ "$report" == "true" ]; then
  echo "----------"
  echo "Progress Step 2: Translation to BED format for CrossMap."
  cat $path_file_temporary_map_source | head -5
  echo "----------"
fi



##########
# 3. Translate genomic coordinates between assemblies of the human genome in
#     CrossMap. Introduce header to table for clarity.
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
$path_file_temporary_map_source \
$path_file_temporary_map_product
# Deactivate Virtual Environment.
deactivate
which python3
# Introduce the same header from the source file to the product file.
# CrossMap does not transfer the original header.
cat $path_file_temporary_map_source | awk 'BEGIN { FS="\t"; OFS="\t" } NR == 1' > $path_file_temporary_map_product_header
cat $path_file_temporary_map_product >> $path_file_temporary_map_product_header

# Report.
if [ "$verbosity" == "true" ] && [ "$report" == "true" ]; then
  echo "----------"
  echo "Progress Step 3: Map between genomic assemblies in CrossMap."
  cat $path_file_temporary_map_product_header | head -5
  echo "----------"
fi



##########
# 4. Organize the table of final genomic coordinates.
echo "IDENTIFIER_MERGE chrom chromStart chromEnd" > $path_file_temporary_merge_product
cat $path_file_temporary_map_product_header | awk 'BEGIN {FS = "\t"; OFS = " "} NR > 1 {
  print $4, $1, $2, $3
}' >> $path_file_temporary_merge_product

# Report.
if [ "$verbosity" == "true" ] && [ "$report" == "true" ]; then
  echo "----------"
  echo "Progress Step 4: Translation of format for merge."
  cat $path_file_temporary_merge_product | head -5
  echo "----------"
fi


##########
# 5. Use common sequential identifier to merge final genomic coordinates to
#     information from the source table.
# Merge text tables by a common identifier in "awk".
# Bash command "join" might also be capable of this type of merge or join, but
# it might require the same identifiers in sort order.
# For computational efficiency, the first file ought to be the subset of the
# second file.
# Table 1 (subset): 4 total columns with merge identifier in column 1.
# Table 2 (complete): 15 total columns with merge identifier in column 1.
# Delimiter: Space
# It is not necessary to print the header row separately.
#echo "IDENTIFIER_MERGE SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT chrom chromStart chromEnd" > $path_file_temporary_product_after_merge
awk 'FNR==NR{a[$1]=$2FS$3FS$4; next} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, a[$1]
}' $path_file_temporary_merge_product $path_file_temporary_merge_source > $path_file_temporary_product_after_merge

cat $path_file_temporary_product_after_merge | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_file_temporary_product_filter
cat $path_file_temporary_product_after_merge | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 18)
    # Skip any rows with incorrect count of column fields.
    next
  else if (($1 == "") || ($1 == "NA"))
    # Skip any rows for which an empty or missing identifier indicates that
    # record in second table did not match first table.
    next
  else
    # Keep record from merge of first table and second table.
    print $0
}' >> $path_file_temporary_product_filter

  # Report.
  if [ "$verbosity" == "true" ] && [ "$report" == "true" ]; then
    echo "----------"
    echo "Progress Step 5: Merge."
    cat $path_file_temporary_product_filter | head -5
    echo "----------"
  fi



##########
# 6. Organize format of product table of GWAS summary statistics.
echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_product_format
cat $path_file_temporary_product_filter | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
  print $2, $16, $18, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15
}' >> $path_file_temporary_product_format

# Compress file format.
gzip -cvf $path_file_temporary_product_format > $path_file_product

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Script:"
  echo $0 # Print full file path to script.
  echo "Translate genomic coordinates in CrossMap between assemblies of the "
  echo "human genome."
  echo "----------"
  echo "path to source file: " $path_file_source
  echo "path to chain file: " $path_file_chain
  echo "path to product file: " $path_file_product
  echo "----------"
  echo "Source Table before translation:"
  echo "Count lines: "
  zcat $path_file_source | wc -l
  zcat $path_file_source | head -5
  echo "----------"
  echo "Product Table after translation:"
  echo "Count lines: "
  zcat $path_file_product | wc -l
  zcat $path_file_product | head -5
  echo "----------"
  echo "----------"
  echo "----------"
fi

# Remove temporary, intermediate files.
rm -r $path_directory_product_temporary



#
