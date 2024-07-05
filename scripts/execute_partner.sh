#!/bin/bash

#chmod u+x script.sh

# TODO: TCW; 5 July 2024
# TODO: This is a very old script that is obsolete.


path_project="/media/tcameronwaller/primary/data/local/work/affiliation/general/partner"
#path_project="/media/tcameronwaller/primary/data/local/work/affiliation/mayo_clinic/projects/bipolar_metabolism"
path_repository="${path_project}/repository"
path_package="${path_repository}/package"
path_python_library="/home/tcameronwaller/project_local/python_library"
path_dock="/home/tcameronwaller/project_local/dock"

# Organize paths to custom package installations.
PYTHONPATH=$path_python_library:$PYTHONPATH
export PYTHONPATH

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

# Collect and organize heritability estimations for metabolites from GWAS
# summary statistics of multiple studies on the human metabolome.
python3 $path_package/interface.py main --path_dock $path_dock --format_tables
