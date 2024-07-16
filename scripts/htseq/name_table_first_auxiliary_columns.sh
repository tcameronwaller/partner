#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 16 July 2024
# Date, last execution or modification: 16 July 2024
# Review: TCW; 16 July 2024
###############################################################################
# Note



###############################################################################
# Organize arguments.

path_file_table_source=${1} # full path to source file for original table
path_file_table_product=${2} # full path to product file for novel table
report=${3} # whether to print reports to terminal

###############################################################################
# Organize paths.

# Initialize file.
rm $path_file_table_product

###############################################################################
# Organize parameters.

###############################################################################
# Execute procedure.

# Extract the original header names of all columns after those that are empty.
# Print the novel header names of the columns that had been empty.
# Print the original header names for the remaining columns.
cat $path_file_table_source | awk -v start=7 'BEGIN {
  FS="\t"; OFS="\t"
} NR==1 {
  columns=""; for (i=start; i<=NF; i++) columns=columns$i OFS; print "identifier_gene", "gene_id", "gene_name", "exon_number", "gene_type", "chromosome", columns
}' > $path_file_table_product
cat $path_file_table_source | awk 'BEGIN {
  FS="\t"; OFS="\t"
} NR > 1 {
  print $0
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



###############################################################################
# End.
