#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, initialization: 17 July 2024
# Date, revision or review: 11 February 2026
###############################################################################
# Note

# This Bash script restores or updates the version of the "partner" repository,
# including its parameters, scripts, package, and subpackages for local
# execution, meaning execution on a local machine rather than on a remote
# server. This restoration does not pull code from public repositories on
# GitHub; rather, it copies code and parameters from a private storage location
# on the local machine to a temporary location for convenient execution.

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
# https://stackoverflow.com/questions/55740123/pip-virtualenv-reset-the-path
# -after-reactivating

###############################################################################
# Organize arguments.

report=${1} # whether to print reports to terminal

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_process=$(<"$path_directory_paths/path_directory_process_local.txt")
path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/data" # restore script does not modify "in_data" for efficiency
path_directory_demonstration="$path_directory_dock/demonstration"
path_directory_parameters="$path_directory_dock/parameters"
path_directory_parameters_private=$(<"$path_directory_paths/path_directory_parameters_private.txt")
path_directory_scripts="$path_directory_process/scripts"
path_directory_package="$path_directory_process/package"
path_directory_repository_partner=$(<"$path_directory_paths/path_directory_repository_partner.txt")

# Initialize directories.
if [ -d $path_directory_demonstration ] || [ -d $path_directory_parameters ] || [ -d $path_directory_scripts ] || [ -d $path_directory_package ] ; then
  # Remove previous versions of code from temporary location for execution.
  rm -rf $path_directory_demonstration
  rm -rf $path_directory_parameters
  rm -rf $path_directory_scripts
  rm -rf $path_directory_package
fi

if [ ! -d $path_directory_process ] || [ ! -d $path_directory_dock ] || [ ! -d $path_directory_demonstration ] || [ ! -d $path_directory_parameters ] || [ ! -d $path_directory_scripts ] || [ ! -d $path_directory_package ] ; then
  # Directory or directories do not already exist.
  # Create directories.
  mkdir -p $path_directory_process
  mkdir -p $path_directory_dock
  mkdir -p $path_directory_demonstration
  mkdir -p $path_directory_parameters
  mkdir -p $path_directory_scripts
  mkdir -p $path_directory_package
fi

################################################################################
# Execute procedure.

##########
# Organize parameters.
cp -r "$path_directory_repository_partner/demonstration" "$path_directory_demonstration/demonstration"
mv "$path_directory_demonstration/demonstration" "$path_directory_demonstration/partner"
cp -r "$path_directory_repository_partner/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/partner"
cp -r $path_directory_parameters_private "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/privacy_security"

##########
# Organize scripts.
cp -r "$path_directory_repository_partner/scripts" "$path_directory_scripts/scripts"
mv "$path_directory_scripts/scripts" "$path_directory_scripts/partner"

##########
# Organize Python packages.
# Package: "partner"
# Hierarchy: child, lower-level subpackage
# Scripts remain within original repository's structure.
# Python code transfers to a subpackage child directory within the parent
# directory of the main package.
cp -r "$path_directory_repository_partner/package" "$path_directory_package/package"
mv "$path_directory_package/package" "$path_directory_package/partner"

##########
# Initialize directory permission.
chmod -R 0777 $path_directory_process

###############################################################################
# Report.
if [ "$report" == "true" ]; then
  echo "----------"
  echo "project: partner"
  echo "script: restore_repository_local.sh"
  echo $0 # Print full file path to script.
  echo "done"
  echo "----------"
fi

###############################################################################
# End.
