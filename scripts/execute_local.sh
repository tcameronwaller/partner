#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

################################################################################
# Activate Python Virtual Environment.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_tools="/media/tcameronwaller/primary/data/local/work/affiliation/general/tools"
path_environment_main="${path_tools}/python/environments/main"
source "${path_environment_main}/bin/activate"
echo "----------"
echo "confirm activation of Python Virtual Environment..."
which python3
sleep 5s

################################################################################

path_process="/media/tcameronwaller/primary/data/local/work/affiliation/general/process_local"
path_dock="$path_process/dock"
path_package="${path_process}/local/local"

# Echo each command to console.
set -x

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_dock ]
then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_dock
fi

# Execute procedure(s).

#python3 $path_package/interface.py main --path_dock $path_dock --scratch

################################################################################
# Deactivate Python Virtual Environment.
deactivate
echo "----------"
echo "confirm deactivation of Python Virtual Environment..."
which python3
