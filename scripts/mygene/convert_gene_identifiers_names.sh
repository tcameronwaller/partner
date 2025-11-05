#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 31 October 2024
# Date, last execution or modification: 24 September 2025
# Review: TCW; 24 September 2025
###############################################################################
# Note

# tool: MyGene.info
# site, home: https://www.mygene.info/
# site, live API: https://www.mygene.info/v3/api
# site, documentation: https://docs.mygene.info/en/latest/index.html
# site, host PyPi: https://pypi.org/project/mygene/

###############################################################################
# Organize arguments.

path_file_source=${1} # full path to source file in text format from which to read identifiers or names of genes
path_file_product_list=${2} # full path to product file to which to write as a list in text format the identifiers or names of genes
path_file_product_table=${3} # full path to product file to which to write as a table in text format the identifiers or names of genes
delimiter_source=${4} # text delimiter between items in source file, not white space
delimiter_product=${5} # text delimiter between items in product file, not white space
type_source=${6} # type of identifiers or names of genes in source query
type_product=${7} # type of identifiers or names of genes in product delivery
species=${8} # name of species, either "human", "mouse", or another relevant option
report=${9} # whether to print reports to terminal

###############################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"

# Files.

# Scripts.
path_file_script_source="${path_directory_scripts}/partner/python/entities_identifiers/mygene_convert_gene_identifiers_names.py"
path_file_script_product="${path_directory_package}/mygene_convert_gene_identifiers_names.py"

# Copy Python script to package directory.
cp $path_file_script_source $path_file_script_product

# Executable handles.
path_environment_main="$path_directory_tools/python/environments/main"
echo $path_environment_main

# Initialize directory.

# Initialize file.
#rm $path_file_product

###############################################################################
# Organize parameters.

# Parameters.
threads=6
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error

###############################################################################
# Activate Python virtual environment.

# Activate Python virtual environment.
source "${path_environment_main}/bin/activate"
# Set paths for local packages and modules.
#echo "Python path variable before update"
#echo $PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$path_directory_package
export PYTHONPATH=$PYTHONPATH:$path_directory_package_partner
export PYTHONPATH=$PYTHONPATH:$path_directory_package_project_main
# Regulate concurrent or parallel process threads on node cores.
# Force Python program (especially SciPy) not to use all available cores on a
# cluster computation node.
export MKL_NUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export OMP_NUM_THREADS=$threads
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "Python virtual environment: main"
  echo "path to Python installation:"
  which python3
  echo "Python path variable:"
  echo $PYTHONPATH
  sleep 1s
  echo "----------"
fi



###############################################################################
# Execute procedure.

# Execute program process in Python.
python3 $path_file_script_product \
$path_file_source \
$path_file_product_list \
$path_file_product_table \
$delimiter_source \
$delimiter_product \
$type_source \
$type_product \
$species \
$report

###############################################################################
# Deactivate Python virtual environment.

# Deactivate Python virtual environment.
deactivate
#which python3

###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: convert_gene_identifiers_names.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
  echo "Convert identifiers or names of genes by query to MyGene.info"
  echo "path to source file: " $path_file_source
  echo "path to product file list: " $path_file_product_list
  echo "path to product file table: " $path_file_product_table
  echo "delimiter in source file: " $delimiter_source
  echo "delimiter in product file: " $delimiter_product
  echo "type of identifiers or names in source: " $type_source
  echo "type of identifiers or names in product: " $type_product
  echo "species: " $species
  echo "----------"
  echo "----------"
  echo "----------"
fi

###############################################################################
# End.
