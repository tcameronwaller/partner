#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, initialization: 4 March 2026
# Date, revision or review: 4 March 2026
###############################################################################
# Note

# Use this script to create a table of parameters for many different features.
# As source, this script reads from file a list of names of features in text
# format with new lines as delimiters.
# The table of parameters is in format to pass to the Python script in file
# "interact_feature_products.py", which controls the calculation of products
# between pairs of features for use as terms for interaction in regression.

###############################################################################
# Organize arguments.

path_file_source_list_features=${1} # full path to file for a list of features
path_file_product_table_parameters=${2} # full path to file for a table of parameters for interaction products
category=${3}
name_interaction_prefix=${4}
feature_first=${5}
features_second_extra=${6}
features_continuity_scale_extra=${7}
date_review=${8}
note=${9}
report=${10} # whether to print reports to terminal

###############################################################################
# Organize parameters.

execution="1"

################################################################################
# Organize paths.

# Simplify name of variable for path.
path_file_product="$path_file_product_table_parameters"

# Initialize file.
rm $path_file_product # caution

###############################################################################
# Execute procedure.

# Write column header names in table to file.
printf "execution\tsequence\tcategory\tname_interaction\t" > $path_file_product
printf "feature_first\tfeature_second\t" >> $path_file_product
printf "features_continuity_scale\t" >> $path_file_product
printf "date_review\tnote\n" >> $path_file_product

# Initialize the sequence counter.
sequence=1

# Collect features from the extra list.
if [[ "$features_second_extra" != "none" ]]; then
  IFS=',' read -r -a array_items <<< "$features_second_extra"
  unset IFS
  # Iterate on features from extra list.
  # Assemble parameters for features from extra list.
  for feature_second in "${array_items[@]}"; do
    # Determine name of feature for interaction product.
    name_interaction="${name_interaction_prefix}_${feature_first}_${feature_second}"
    # Write rows in table to file.
    printf "${execution}\t${sequence}\t${category}\t${name_interaction}\t" >> $path_file_product
    printf "${feature_first}\t${feature_second}\t" >> $path_file_product
    printf "${features_continuity_scale_extra}\t" >> $path_file_product
    printf "${date_review}\t${note}\n" >> $path_file_product
    # Increment sequence.
    ((sequence++))
  done
fi

# Read and collect features from the main list in the source file.
# Initialize array.
array_items=()
# Determine whether file exists for list of features.
path_file_source=$path_file_source_list_features
if [[ -f $path_file_source ]]; then
  # Read text items from file in text format with newline delimiters.
  #readarray -t array_items < $path_file_source
  #IFS='\n' read -r -a array_items <<< $path_file_source
  while IFS=$'\n' read -r -a item
  do
    # Determine whether the item is an empty string.
    if [[ -n "$item" ]]; then
      # Collect.
      array_items+=("${item}")
    fi
  done < <(tail -n +0 "${path_file_source}"; echo) # append new line to tail to ensure read of last line
  unset IFS
  # Iterate on features from main list.
  # Assemble parameters for features from main list.
  for feature_second in "${array_items[@]}"; do
    # Determine name of feature for interaction product.
    name_interaction="${name_interaction_prefix}_${feature_first}_${feature_second}"
    # Write rows in table to file.
    printf "${execution}\t${sequence}\t${category}\t${name_interaction}\t" >> $path_file_product
    printf "${feature_first}\t${feature_second}\t" >> $path_file_product
    printf "${feature_second}\t" >> $path_file_product
    printf "${date_review}\t${note}\n" >> $path_file_product
    # Increment sequence.
    ((sequence++))
  done
  unset IFS
fi

# Report.
if [[ "$report" == "true" ]]; then
  # Organize information for report.
  count_items=${#array_items[@]}
  index_array_maximum=$(($count_items - 1))
  # Report information.
  echo "----------"
  echo "count of items: " $count_items
  echo "first item: " ${array_items[0]} # notice base-zero indexing
  echo "last item: " ${array_items[$index_array_maximum]}
  echo "----------"
fi

###############################################################################
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: assemble_table_parameters_interaction.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
  echo "----------"
  echo "----------"
fi



###############################################################################
# End.
