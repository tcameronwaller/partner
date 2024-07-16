#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 15 July 2024
# Date, last execution or modification: 15 July 2024
# Review: TCW; 15 July 2024
###############################################################################
# Note

# TODO: TCW; 16 July 2024
# TODO: This script is not quite ready. I still need to adjust the count of
# columns in each novel supplement table that needs to be merge to main.


###############################################################################
# Organize arguments.

path_file_table_main=${1} # full path to source file for original main table to which to merge columns from novel supplement table
path_file_table_supplement=${2} # full path to source file for novel supplement table to merge with original main table
path_file_table_merge=${3} # full path to product file for table after merge
report=${4} # whether to print reports to terminal

###############################################################################
# Organize paths.

# Remove any previous version of the product file.
rm $path_file_table_merge

###############################################################################
# Organize parameters.

###############################################################################
# Execute procedure.

##########
# Merge novel supplement table to original main table.
# Novel supplement table must have specific count of columns, all of which
# transfer to the original main table.
# Original main table can have any count of columns, all of which remain.
# Merge text tables by a common identifier in "awk".
# Bash command "join" might also be capable of this type of merge or join, but
# it might require the same identifiers in sort order.
# These merges treat Table 2 as the priority table.
# The merge of Table 1 with Table 2 will include all records from Table 2 but
# only those records from Table 1 that have matching identifiers.
# Any identifiers in Table 2 that are not also in Table 1 will correspond to
# missing values across all columns from Table 1.
# These merges print each row of Table 2 regardless of whether there is a record
# with matching identifier in Table 1.

# Process explanation:
# 1. GNU awk reads lines from the first file (FNR==NR) into an array in memory
# (array "a") that uses values from the first column of the first file (column
# "identifier_merge") as a hashable index (a[$1]=...).
# 2. GNU awk then proceeds to the second file and prints line by line while
# including any of the hashable lines from the first file that have a matching
# index, this time from the first column of the second file (column
# "identifier_merge").

# Reference:
# https://stackoverflow.com/questions/32481877/what-are-nr-and-fnr-and-what-does-nr-fnr-imply



##########
# Merge
# Table 1: novel supplement table
# Table 1: 7 total columns with merge identifiers in column 1
# Table 2: original main table
# Table 2: __ to many columns with merge identifiers in column 1
# Merged table: many total columns with merge identifiers in column 1 and preceding block of merged columns from table 1
# Delimiter: Tab ("\t")
# It is not necessary to print the header row separately as long as the column of identifiers has the same header name in both tables.
# Note that AWK interprets a single space delimiter (FS=" ") as any white space.
#echo "identifier_merge custom column names" > $path_file_table_merge
awk 'BEGIN {
  FS="\t"; OFS="\t"
} FNR==NR {
  a[$1]=$2FS$3FS$4FS$5FS$6FS$7; next
} {
  if(a[$1]==""){a[$1]="NA"FS"NA"FS"NA"FS"NA"FS"NA"FS"NA"}; print $0, a[$1]
} END {
  delete a
}' $path_file_table_supplement $path_file_table_main > $path_file_table_merge



##########
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Source main table:"
  head -10 $path_file_table_main
  echo "----------"
  echo "Source supplement table:"
  head -10 $path_file_table_supplement
  echo "----------"
  echo "Product table after merge:"
  head -10 $path_file_table_merge
  echo "----------"
  echo "----------"
  echo "----------"
fi



###############################################################################
# End.
