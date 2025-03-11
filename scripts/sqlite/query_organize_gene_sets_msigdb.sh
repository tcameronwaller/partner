#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 7 March 2025
# Date, last execution or modification: 7 March 2025
# Review: TCW; 7 March 2025
################################################################################
# Note

# Access the Molecular Signatures Database (MSigDB) as an SQLite database.
# site: https://www.gsea-msigdb.org/gsea/downloads.jsp
# wget https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2024.1.Hs/msigdb_v2024.1.Hs.db.zip
# unzip msigdb_v2024.1.Hs.db.zip -d msigdb_v2024.1.Hs.db
# mv msigdb_v2024.1.Hs.db msigdb_sqlite
# cp -r msigdb_v2024.1.Hs.db "$path_directory_dock/in_data/msigdb"

################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"

path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_parameters_private="$path_directory_dock/in_parameters_private"

path_directory_product="$path_directory_dock/out_msigdb_sqlite"

# Files.
path_file_source_msigdb="${path_directory_data}/msigdb/msigdb_sqlite/msigdb_v2024.1.Hs.db"
path_file_table_selection="${path_directory_parameters_private}/age_exercise/reference/table_msigdb_gene_sets_selection.tsv"
# "species" "collection_name" "standard_name"

path_file_product="${path_directory_dock}/out_msigdb_sqlite/msigdb_selection.txt"

# Scripts.
path_file_script_source="${path_directory_scripts}/partner/python/script_sqlite_msigdb_query_organize_gene_sets.py"
path_file_script_product="${path_directory_package}/script_sqlite_msigdb_query_organize_gene_sets.py"
# Copy Python script to package directory.
cp $path_file_script_source $path_file_script_product
# Executable handles.
path_environment_main="$path_directory_tools/python/environments/main"
echo $path_environment_main

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



###############################################################################
# Organize parameters.

# Parameters.
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error
threads=6
report="true"
identifier_type="symbol_hugo" # "symbol_hugo", "entrezgene", "ensembl.gene",
delimiter_name="tab" # "newline", "tab", "space", "semicolon", "colon", "comma",

###############################################################################
# Activate Python virtual environment.

# Activate Python virtual environment.
source "${path_environment_main}/bin/activate"

# Set paths for local packages and modules.
export OLD_PYTHONPATH="$PYTHONPATH"
export PYTHONPATH="$PYTHONPATH:$path_directory_package"

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



################################################################################
# Execute procedure.

# Execute program process in Python.
python3 $path_file_script_product \
$path_file_source_msigdb \
$path_file_table_selection \
$path_file_product \
$identifier_type \
$delimiter_name



###############################################################################
# Deactivate Python virtual environment.

# Restore paths.
export PYTHONPATH="$OLD_PYTHONPATH"

# Deactivate Python virtual environment.
deactivate
#which python3



###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "script: query_organize_gene_sets_msigdb.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

##########
# Remove directory of temporary, intermediate files.
#rm -r $path_directory_temporary



###############################################################################
# End.



#
