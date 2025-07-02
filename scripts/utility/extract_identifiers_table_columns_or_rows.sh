#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 2 July 2025
# Date, last execution or modification: 2 July 2025
# Review: 2 July 2025
###############################################################################
# Note

# TODO: TCW; 2 July 2025

# This script is a place holder and a plan.
# The script will accept arguments...
# path to directory and file for table in tab-delimited text format
# whether to extract identifiers from columns or rows
# index of column or row from which to extract identifiers
# delimiter for the extracted identifiers
# path to directory and file for the product list in text format


###############################################################################
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: extract_identifiers_table_columns_or_rows.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
  #echo "path to file for source table: " $path_file_table_source
  #echo "path to file for product table: " $path_file_table_product
  echo "----------"
  echo "----------"
  echo "----------"
fi

###############################################################################
# End.
