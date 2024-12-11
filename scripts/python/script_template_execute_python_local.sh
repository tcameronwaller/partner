#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: ___
# Date, last execution or modification: ___
# Review: ___
###############################################################################
# Note



###############################################################################
# Organize arguments.



###############################################################################
# Organize paths.

# Project.
project_main="name_of_repository"
interface_subparser="name_of_subparser_in_interface"

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"
path_directory_package_partner="$path_directory_package/partner"
path_directory_package_project_main="$path_directory_package/${project_main}"

path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_parameters="$path_directory_dock/in_parameters/${project_main}"
path_directory_parameters_private="$path_directory_dock/in_parameters_private/${project_main}"
#path_directory_source="${path_directory_dock}/${project_main}/source"
#path_directory_product="${path_directory_dock}/${project_main}/product"
#stamp_date=$(date +%Y-%m-%d)
#path_directory_temporary="${path_directory_product}/temporary_${stamp_date}" # hopefully unique

# Files.

# Scripts.

# Executable handles.
path_environment_main="$path_directory_tools/python/environments/main"
echo $path_environment_main
# Initialize directory.
#rm -r $path_directory_product # caution
#mkdir -p $path_directory_product
#mkdir -p $path_directory_temporary

# Initialize file.



###############################################################################
# Organize parameters.

# Parameters.
threads=6
report="true"
#set -x # enable print commands to standard error
set +x # disable print commands to standard error
#set -v # enable print input to standard error
set +v # disable print input to standard error



###############################################################################
# Activate Python virtual environment.

# Activate Python virtual environment.
source "${path_environment_main}/bin/activate"

# Set paths for local packages and modules.
# Normally, activation of the Python virtual environment will update the
# "PYTHONHOME" and "PYTHONPATH" variables as necessary to point to the proper
# installations of Python and its packages.
# Normally it is only necessary to update the "PYTHONPATH" variable to point to
# custom packages and modules in a special location on the system.

#export OLD_PYTHONHOME="$PYTHONHOME"
export OLD_PYTHONPATH="$PYTHONPATH"
#export PYTHONHOME="$PYTHONHOME:$VIRTUAL_ENV"
#export PYTHONPATH="$PYTHONPATH:$VIRTUAL_ENV"
#export PYTHONHOME="$PYTHONHOME:${path_environment_main}"
#export PYTHONPATH="$PYTHONPATH:${path_environment_main}"
#export PYTHONHOME="$PYTHONHOME:$path_directory_package"
#export PYTHONHOME="$PYTHONHOME:$VIRTUAL_ENV/lib/python-3.13/site-packages"
#export PYTHONPATH="$PYTHONPATH:$VIRTUAL_ENV/lib/python-3.13/site-packages"
#export PYTHONHOME="$PYTHONHOME:${path_environment_main}/lib/python3.13/site-packages"
#export PYTHONPATH="$PYTHONPATH:${path_environment_main}/lib/python3.13/site-packages"
export PYTHONPATH="$PYTHONPATH:$path_directory_package"
#export PYTHONPATH=$PYTHONPATH:$path_directory_package_partner
#export PYTHONPATH=$PYTHONPATH:$path_directory_package_project_main

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
  echo "VIRTUAL_ENV variable:"
  echo $VIRTUAL_ENV
  echo "PYTHONHOME variable:"
  echo $PYTHONHOME
  echo "PYTHONPATH variable:"
  echo $PYTHONPATH
  sleep 1s
  echo "----------"
fi



###############################################################################
# Execute procedure.

# Execute program process in Python.
python3 $path_directory_package_project_main/interface.py \
$interface_subparser \
--module_procedure \
--path_directory_dock $path_directory_dock

###############################################################################
# Deactivate Python virtual environment.

# Restore paths.
#export PYTHONHOME="$OLD_PYTHONHOME"
export PYTHONPATH="$OLD_PYTHONPATH"

# Deactivate Python virtual environment.
deactivate
#which python3



###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "project: ${project_main}"
  echo "routine:: ${interface_subparser}"
  echo "procedure: 3_differential_expression"
  echo "script: 3_interaction.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

##########
# Remove directory of temporary, intermediate files.
#rm -r $path_directory_temporary



###############################################################################
# End.
