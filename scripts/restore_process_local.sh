#!/bin/bash

#chmod u+x script.sh

###########################################################################
# Organize script parameters.

project="local"

path_process="/media/tcameronwaller/primary/data/local/work/affiliation/general/process_local"
path_repository="$path_process/local"
path_parameters="$path_process/dock/parameters"

path_partner_source="/media/tcameronwaller/primary/data/local/work/affiliation/general/partner/repository"
path_partner="$path_process/partner"

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous versions of the repositories..."
rm -rf $path_partner
rm -r $path_parameters

##########
# Organize and restore supplemental sub-repositories.

# Repository: partner
# Scripts remain within original repository's structure.
# Python code transfers to sub-package.
cp -r "$path_partner_source" "$path_partner"
mv "$path_partner/package" "$path_partner/partner"
cp -r "$path_partner/partner" "$path_repository/${project}/partner"

##########
# Organize and restore parameters.

mkdir -p $path_parameters
cp -r "$path_repository/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/${project}"
cp -r "$path_partner/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/partner"
