#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 26 June 2024
# Date, last execution or modification: 26 June 2024
# Review: TCW; 26 June 2024
################################################################################
# Note

# This Bash script pulls code from GitHub to restore (update) the version of the
# 'partner' project repository, which includes Bash scripts, Python scripts, and
# a Python package. Since the 'partner' Python package is designed to execute as
# a subpackage under the management of a higher level package, this script
# sets up a dummy top-level package from which to execute Python scripts that in
# turn use functionality from the 'partner' Python subpackage.


################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_process=$(<"./paths/psychiatry/process_aai.txt")
path_directory_dock="$path_directory_process/dock"
path_directory_data="$path_directory_dock/in_data"
path_directory_parameters="$path_directory_dock/in_parameters"
path_directory_package="$path_directory_process/package_top_dummy"

path_directory_partner_product="$path_directory_process/partner"

# Initialize directories.

if [ -d $path_directory_parameters ] || [ -d $path_directory_package ] || [ -d $path_directory_partner_product ] ; then
  # Remove previous versions of code from temporary location for execution.
  rm -rf $path_directory_parameters
  rm -rf $path_directory_package
  rm -rf $path_directory_partner_product
fi

if [ ! -d $path_directory_process ] || [ ! -d $path_directory_dock ] || [ ! -d $path_directory_parameters ] || [ ! -d $path_directory_package ] ; then
  # Directory or directories do not already exist.
  # Create directories.
  mkdir -p $path_directory_process
  mkdir -p $path_directory_package
  mkdir -p $path_directory_dock
  mkdir -p $path_directory_parameters
fi

################################################################################
# Execute procedure.

# Echo each command to console.
set -x

##########
# Access and organize current version of repository.
# Repository: "partner"
# Hierarchy: child, lower-level subpackage
# Scripts remain within original repository's structure.
# Python code transfers to a subpackage child directory within the parent
# directory of the main package.
cd $path_directory_process
wget https://github.com/tcameronwaller/partner/archive/main.zip
unzip main.zip
rm main.zip
mv partner-main $path_directory_partner_product
cp -r "$path_directory_partner_product/package" "$path_directory_package"
mv "$path_directory_package/package" "$path_directory_package/partner"

##########
# Copy and organize current version of parameters.
cp -r "$path_directory_partner_product/parameters" "$path_directory_parameters/parameters"
mv "$path_directory_parameters/parameters" "$path_directory_parameters/partner"

##########
# Create dummy top-level package
cp "$path_directory_package/partner/__init__.py" "$path_directory_package/__init__.py"

##########
# Initialize directory permission.
chmod -R 0777 $path_directory_process

################################################################################
# Report.
echo "----------"
echo "Script complete:"
echo $0 # Print full file path to script.
echo "restore_repository_partner_dummy_top.sh"
echo "----------"

################################################################################
# End.
