#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 16 July 2024
# Date, last execution or modification: 17 July 2024
# Review: TCW; 17 July 2024
###############################################################################
# Note



###############################################################################
# Organize arguments.

path_file_table_source=${1} # full path to source file for original table
path_file_table_product=${2} # full path to product file for novel table
report=${3} # whether to print reports to terminal

###############################################################################
# Organize paths.

# Directories.
stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_table_product .tsv)"
path_directory_product="$(dirname $path_file_table_product)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique

# Files.
path_file_temporary_1="${path_directory_temporary}/${name_base_file_product}_temporary_1.tsv"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_temporary

# Initialize file.
rm $path_file_temporary_1
rm $path_file_table_product

###############################################################################
# Organize parameters.

###############################################################################
# Execute procedure.

# Copy table to product file.
cp $path_file_table_source $path_file_temporary_1

# Replace problematic delimiters.
# Rather than a true tab delimiter, the tables in text format have delimiters
# of literal "\t" string characters between columns in each row.
sed -i 's%\\t%;%g' $path_file_temporary_1

# Extract the original header names of all columns after those that are empty.
# Print the novel header names of the columns that had been empty.
# Print the original header names for the remaining columns.
cat $path_file_temporary_1 | awk -v start=7 'BEGIN {
  FS=";"; OFS="\t"
} NR==1 {
  columns=""; for (i=start; i<=NF; i++) columns=columns$i OFS; print "identifier_gene", "gene_id", "gene_name", "exon_number", "gene_type", "chromosome", columns
}' > $path_file_table_product
# To change the delimiters, print each column within each row.
cat $path_file_temporary_1 | awk -v start=1 'BEGIN {
  FS=";"; OFS="\t"
} NR > 1 {
  columns=""; for (i=start; i<=NF; i++) columns=columns$i OFS; print columns
}' >> $path_file_table_product



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Source original table:"
  head -10 $path_file_table_source
  echo "----------"
  echo "Product novel table:"
  head -10 $path_file_table_product
  echo "----------"
  echo "----------"
  echo "----------"
fi

##########
# Remove directory of temporary, intermediate files.
rm -r $path_directory_temporary

###############################################################################
# End.
