#!/bin/bash

#chmod u+x script.sh

# Read in private file paths.
echo "read private file path variables..."
cd ~/path
path_temporary=$(<"./process_temporary.txt")
path_project=$(<"./project_sexy_alcohol.txt")
path_repository=$(<"./sexy_alcohol.txt")

# Echo each command to console.
set -x

# download python
# https://www.python.org/downloads/
# /usr/bin/python3

# Install pip unless already installed with python.
# https://pip.pypa.io/en/stable/installing/
#curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
#python get-pip.py --target <directory>
#pip3 install --upgrade pip

# Install SciPy.
pip3 install scipy --upgrade --target <path>

# Install NumPy.
pip3 install numpy --upgrade --target <path>

# Install Pandas.
pip3 install pandas --upgrade --target <path>

# Install SkLearn.
pip3 install sklearn --upgrade --target <path>

# Install StatsModels.
pip3 install statsmodels --upgrade --target <path>

# Install StatsModels.
pip3 install matplotlib --upgrade --target <path>


# tutorial on paths
#https://www.devdungeon.com/content/python-import-syspath-and-pythonpath-tutorial#toc-5
