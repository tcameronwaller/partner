#!/bin/bash

#chmod u+x script.sh

path_python="/usr/bin/python3"
path_python_library="/home/tcameronwaller/project_local/python_library"

# Echo each command to console.
set -x

cd ~

# Python.
python3 --version
#wget "https://www.python.org/ftp/python/3.9.4/Python-3.9.4.tgz"
#sudo apt update
#sudo apt install python3.9

# Pip.
#python3 -m pip --version
#python3 -m ensurepip --upgrade
#sudo apt install python3-pip
#python3 -m pip3 install --upgrade pip3
pip3 --version

pip3 install --target=${path_python_library} --upgrade numpy
pip3 install --target=${path_python_library} --upgrade scipy
pip3 install --target=${path_python_library} --upgrade testresources
pip3 install --target=${path_python_library} --upgrade pandas
pip3 install --target=${path_python_library} --upgrade sklearn
pip3 install --target=${path_python_library} --upgrade statsmodels
pip3 install --target=${path_python_library} --upgrade matplotlib
pip3 install --target=${path_python_library} --upgrade networkx
