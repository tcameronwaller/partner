#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 4 April 2025
# Date, last execution or modification: 4 April 2025
# Review: 4 April 2025
###############################################################################
# Note

# Use this script to create a table of parameters for many different response
# features but otherwise identical predictor features and other parameters.
# As source, this script reads from file a list of names of response features
# in text format with new lines as delimiters.
# The table of parameters is in format to pass to the Python script in file
# "script_drive_regressions_from_table_parameters.py", which drives multiple
# regressions in parallel.



###############################################################################
# Organize arguments.



################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"
path_directory_package_partner="$path_directory_package/partner"

path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_demonstration="$path_directory_dock/in_demonstration"
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_parameters_private="$path_directory_dock/in_parameters_private"

path_directory_source="${path_directory_demonstration}/partner"
path_directory_product="${path_directory_dock}/out_regression/demonstration"
#stamp_date=$(date +%Y-%m-%d)
#path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Files.
path_file_source="${path_directory_source}/list_regression_responses.txt"
path_file_product="${path_directory_product}/table_regression_parameters_automatic.tsv"

# Initialize directory.
mkdir -p $path_directory_product

# Initialize file.
rm $path_file_product # caution

###############################################################################
# Organize parameters.

# Parameters.
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

# Common parameters for all instances of parameters for regression.
execution="1"
sequence=1
group="group_automatic"
instance="instance_automatic"
selection_observations="group:group_1,group_2,group_3,group_4;sex:female,male"
type_regression="continuous_ols"
formula_text="response ~ predictor_fixed_1 + predictor_fixed_2 + predictor_fixed_3 + sex_y"
#feature_response="${response}" # this parameter varies
features_predictor_fixed="predictor_fixed_1,predictor_fixed_2,predictor_fixed_3,sex_y"
features_predictor_random="none"
features_continuity_scale="predictor_fixed_1,predictor_fixed_2"
identifier_observations="observation"
method_scale="z_score"
data_path_directory="dock,in_data,regression,demonstration"
data_file="table_regression_data.tsv"
review="2025-04-04"
note="a script prepared this table of parameters automatically"



###############################################################################
# Execute procedure.

# Read list or array of response features from file.
##########
# Read text items from file to array.
# Parameters.
delimiter_source="newline" # "newline", "tab", "\n", "\t", ";", ":", ",", not " "
# Initialize array.
items_source=() # lines
# Read text items from file using delimiters such as new line.
input=$path_file_source
if [[ "$delimiter_source" == "newline" ]]; then
  while IFS=$'\n' read -r -a item
  do
  # Report.
  if [ "$report" == "true" ]; then
    echo "----------"
    echo "item: ${item}"
    echo "----------"
  fi
  # Collect.
  items_source+=("${item}")
done < <(tail -n +0 "${input}"; echo) # append new line to tail to ensure read of last line
fi

# Alternative.

# Read items in array.
readarray -t items_source < $path_file_source
count_items=${#items_source[@]}
index_array_maximum=$(($count_items - 1))

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "count of items: " $count_items
  echo "first item: " ${items_source[0]} # notice base-zero indexing
  echo "last item: " ${items_source[$index_array_maximum]}
  echo "----------"
fi

# Write column header names in table to file.
printf "execution\tsequence\tgroup\tinstance\t" > $path_file_product
printf "selection_observations\ttype_regression\t" >> $path_file_product
printf "formula_text\tfeature_response\t" >> $path_file_product
printf "features_predictor_fixed\tfeatures_predictor_random\t" >> $path_file_product
printf "features_continuity_scale\tidentifier_observations\t" >> $path_file_product
printf "method_scale\tdata_path_directory\tdata_file\t" >> $path_file_product
printf "review\tnote\n" >> $path_file_product
# Write rows in table to file.
# Iterate on response features in array.
# For each response feature, create a new row of parameters in the table.
for item_source in "${items_source[@]}"; do
  printf "${execution}\t${sequence}\t${group}\t${instance}\t" >> $path_file_product
  printf "${selection_observations}\t${type_regression}\t" >> $path_file_product
  printf "${formula_text}\t${item_source}\t" >> $path_file_product
  printf "${features_predictor_fixed}\t${features_predictor_random}\t" >> $path_file_product
  printf "${features_continuity_scale}\t${identifier_observations}\t" >> $path_file_product
  printf "${method_scale}\t${data_path_directory}\t${data_file}\t" >> $path_file_product
  printf "${review}\t${note}\n" >> $path_file_product
  # Increment sequence.
  ((sequence++))
done



###############################################################################
# Report.

if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: template_create_parameter_table.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
  echo "----------"
  echo "----------"
fi



###############################################################################
# End.
