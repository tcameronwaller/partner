#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 17 July 2024
# Date, last execution or modification: 14 October 2024
# Review: TCW; 14 October 2024
###############################################################################
# Note

# This Bash script restores the version of the ${project_main} repository,
# including its parameters, scripts, package, and subpackages for local
# execution, meaning execution on a local machine rather than a remote server.
# This restoration does not pull code from public repositories on GitHub;
# rather, it copies code from a private, more permanent storage location on the
# local machine to a more temporary location for convenient process execution.

# It is unnecessary to copy repositories to the working directory for
# execution; however, it is necessary for the parent directory of each Python
# package to have the appropriate name. In order to store Python code within
# the "package" subdirectories of repositories, it is convenient to copy and
# rename the package directories.

# It is possible for the paths specific to a Python virtual environment to
# become corrupted. It might be necessary to examine the paths in the
# environment's "activate" script. In particular, any change to the names of
# directories in the file system's path to the virtual environment will disrupt
# the environment since the "activate" script uses absolute paths. The files
# for the relevant scripts have names "activate", "activate.csh", and
# "activate.fish". If correcting the paths in these activation scripts is not
# sufficient, then it might be necessary to re-create the environment.
# https://stackoverflow.com/questions/55740123/pip-virtualenv-reset-the-path-after-reactivating

###############################################################################
# Organize arguments.

project_main=${1} # name of main project with repository, package, and parameters
report=${2} # whether to print reports to terminal

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_dock="$path_directory_process/dock"
#path_directory_data="$path_directory_dock/in_data" # restore script does not modify "in_data" for efficiency
path_directory_demonstration="$path_directory_dock/in_demonstration"
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_parameters_private_source=$(<"$path_directory_paths/path_directory_parameters_private_${project_main}.txt")
path_directory_parameters_private="$path_directory_dock/in_parameters_private"
path_directory_repository_partner=$(<"$path_directory_paths/path_directory_repository_partner.txt")
path_directory_repository_project_main=$(<"$path_directory_paths/path_directory_repository_${project_main}.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"
path_directory_package_partner_source="$path_directory_repository_partner/package"
path_directory_package_project_main_source="$path_directory_repository_project_main/package"
path_directory_package_partner_product="$path_directory_package/partner"
path_directory_package_project_main_product="$path_directory_package/${project_main}"

# Initialize directories.
if [ -d $path_directory_demonstration ] || [ -d $path_directory_parameters ] || [ -d $path_directory_parameters_private ] || [ -d $path_directory_scripts ] || [ -d $path_directory_package ] ; then
  # Remove previous versions of code from temporary location for execution.
  rm -rf $path_directory_demonstration
  rm -rf $path_directory_parameters
  rm -rf $path_directory_parameters_private
  rm -rf $path_directory_scripts
  rm -rf $path_directory_package
fi

if [ ! -d $path_directory_process ] || [ ! -d $path_directory_package ] || [ ! -d $path_directory_dock ] || [ -d $path_directory_demonstration ] || [ ! -d $path_directory_parameters ] || [ ! -d $path_directory_parameters_private ] || [ ! -d $path_directory_scripts ] ; then
  # Directory or directories do not already exist.
  # Create directories.
  mkdir -p $path_directory_process
  mkdir -p $path_directory_package
  mkdir -p $path_directory_dock
  mkdir -p $path_directory_demonstration
  mkdir -p $path_directory_parameters
  mkdir -p $path_directory_parameters_private
  mkdir -p $path_directory_scripts
fi

################################################################################
# Execute procedure.

##########
# Organize parameters.
cp -r "$path_directory_repository_partner/demonstration" "$path_directory_demonstration/demonstration"
mv "$path_directory_demonstration/demonstration" "$path_directory_demonstration/partner"
cp -r "$path_directory_repository_partner/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/partner"
cp -r "$path_directory_repository_project_main/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/${project_main}"
cp -r $path_directory_parameters_private_source $path_directory_dock
mv "${path_directory_dock}/parameters" "${path_directory_dock}/${project_main}"
mv "${path_directory_dock}/${project_main}" $path_directory_parameters_private

##########
# Organize scripts.
cp -r "$path_directory_repository_partner/scripts" "$path_directory_scripts/scripts"
mv "$path_directory_scripts/scripts" "$path_directory_scripts/partner"
cp -r "$path_directory_repository_project_main/scripts" "$path_directory_scripts/scripts"
mv "$path_directory_scripts/scripts" "$path_directory_scripts/${project_main}"

##########
# Organize Python packages.
# Package: "partner"
# Hierarchy: child, lower-level subpackage
# Scripts remain within original repository's structure.
# Python code transfers to a subpackage child directory within the parent
# directory of the main package.
if true; then
  cp -r "$path_directory_package_partner_source" "$path_directory_package"
  mv "$path_directory_package/package" "$path_directory_package_partner_product"
fi
# Package: "${project_main}"
# Hierarchy: child, lower-level subpackage
# Scripts remain within original repository's structure.
# Python code transfers to a subpackage child directory within the parent
# directory of the main package.
if true; then
  cp -r "$path_directory_package_project_main_source" "$path_directory_package"
  mv "$path_directory_package/package" "$path_directory_package_project_main_product"
fi

##########
# Initialize directory permission.
chmod -R 0777 $path_directory_process

###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "project: ${project_main}"
  echo "script: restore_repository_local.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

###############################################################################
# End.
