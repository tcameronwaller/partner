#!/bin/bash

#chmod u+x script.sh

###########################################################################
# Organize script parameters.

project="local"

path_process="/media/tcameronwaller/primary/data/local/work/affiliation/general/process_local"
path_repository="$path_process/local"
path_parameters="$path_process/dock/parameters"

path_promiscuity_source="/media/tcameronwaller/primary/data/local/work/affiliation/general/promiscuity/repository"
path_promiscuity="$path_process/promiscuity"

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous versions of the repositories..."
rm -rf $path_promiscuity
rm -r $path_parameters

##########
# Organize and restore supplemental sub-repositories.

# Repository: promiscuity
# Scripts remain within original repository's structure.
# Python code transfers to sub-package.
cp -r "$path_promiscuity_source" "$path_promiscuity"
mv "$path_promiscuity/package" "$path_promiscuity/promiscuity"
cp -r "$path_promiscuity/promiscuity" "$path_repository/${project}/promiscuity"

##########
# Organize and restore parameters.

mkdir -p $path_parameters
cp -r "$path_repository/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/${project}"
cp -r "$path_promiscuity/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/promiscuity"
