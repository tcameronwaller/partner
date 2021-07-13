#!/bin/bash
#chmod u+x script.sh

################################################################################
# Paths

# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_python="${path_tools}/python"
path_python_396="${path_python}/python_3.9.6"
path_python_2718="${path_python}/python_2.7.18"

path_environment_main="${path_tools}/python/environments/main"
path_environment_ldsc="${path_tools}/python/environments/ldsc"


# Initialize directories.
mkdir -p $path_python
mkdir -p $path_python_396
mkdir -p $path_python_2718



################################################################################
# Python installation
# There is a tool, "pyenv" to help manage installations of multiple versions of
# Python. This tool uses "shims" in the Path variable to determine which version
# to use. The same can be accomplished with less over-head by installing each
# version of Python from source to an explicit, unique path. Specify this
# explicit path for the Python version when creating the Virtual Environment.
# Here is a helpful thread on Stack Exchange about management of Python
# versions.
# https://unix.stackexchange.com/questions/190794/uninstall-python-installed-by-compiling-source
# Anthon recommended using "make install" with unique prefix paths instead
# of using "make altinstall". Anthon also recommended directing these
# version-specific Python installations to a directory other than the operating
# system's versions of Python (/usr/bin... or /usr/local/bin...) to avoid
# conflicts.
# Then point each new Virtual Environment to the specific installation version
# to use.
# Here is another helpful description of this installation of multiple Python
# versions to unique paths.
# https://www.benhepworth.com/blog/2016/05/18/install-python-to-separate-directory-on-linux-in-5-easy-steps/
# It is possible to use a shebang to specify the Python installation to use in
# scripts.
# #! $path_python_396

# TODO: It MIGHT be necessary to install dependencies BEFORE installing Python2 and Python3 from source...

# Latest
# Python 3.9.6, release date: 28 June 2021
# https://www.python.org/downloads/release/python-396/
# https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz
cd $path_python_396
wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz # TCW, 6 July 2021
tar -xzvf Python-3.9.6.tgz
#find ~/python -type d | xargs chmod 0755
cd Python-3.9.6
./configure --prefix=$path_python_396 # TCW, 6 July 2021
make # TCW, 6 July 2021
make test # TCW, 6 July 2021
# 409 tests OK.
# 3 tests failed: test_pathlib, test_urllib2net, test_zipfile
make install # TCW, 6 July 2021

"${path_python_396}/bin/python3" -V # "Python 3.9.6", TCW, 6 July 2021

# Remove installation.
#rm -rf $path_python_396

# Python 2.7.18, release date: 20 April, 2020
# https://www.python.org/downloads/release/python-2718/
# https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz
cd $path_python_2718
wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz # TCW, 6 July 2021
tar -xzvf Python-2.7.18.tgz
#find ~/python -type d | xargs chmod 0755
cd Python-2.7.18
./configure --prefix=$path_python_2718 --with-ensurepip=install # TCW, 6 July 2021
make # TCW, 6 July 2021
make test # TCW, 6 July 2021
# 365 tests OK.
# 1 test failed: test_locale
make install # TCW, 6 July 2021

"${path_python_2718}/bin/python2" -V # Python 2.7.18, TCW, 6 July 2021

# Remove installation.
#rm -rf $path_python_2718



################################################################################
# Pip
# Install pip for specific installations of Python.

"${path_python_396}/bin/python3" -m pip install --user --upgrade pip
"${path_python_396}/bin/python3" -m pip --version # pip 21.1.3, TCW, 6 July 2021

"${path_python_2718}/bin/python2" -m pip install --user --upgrade pip
"${path_python_2718}/bin/python2" -m pip --version # pip 20.3.4, TCW, 6 July 2021



################################################################################
# Python virtual environments
# If you activate a Virtual Environment before installing Python packages to a
# specific path, then the Virtual Environment will remember those paths within
# its own "PYTHONPATH" variable.

# Python standard module "venv" is the preferred Virtual Environment manager for
# Python version 3.3 or higher.
# https://docs.python.org/3/tutorial/venv.html
# [path to specific python3 installation] -m venv [path to new virtual environment and its name to create]
# source [virtual environment name]/bin/activate
# deactivate

"${path_python_396}/bin/python3" -m venv $path_environment_main

# Activate virtual environment.
source "${path_environment_main}/bin/activate"
which python3
pip3 install --target=${path_python_library} --upgrade numpy
pip3 install --target=${path_python_library} --upgrade scipy
pip3 install --target=${path_python_library} --upgrade testresources
pip3 install --target=${path_python_library} --upgrade pandas
pip3 install --target=${path_python_library} --upgrade sklearn
pip3 install --target=${path_python_library} --upgrade statsmodels
pip3 install --target=${path_python_library} --upgrade matplotlib
pip3 install --target=${path_python_library} --upgrade networkx
deactivate
which python3

# Delete virtual environment.
# rm -rf [virtual environment path and name] # remove virtual environment



# Python standard module "virualenv" is the preferred Virtual Environment
# manager for previous versions of Python, including Python 2
# https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/
# Install manager:
# [path to python] -m pip install --user --upgrade virtualenv
# Create new Virtual Environment:
# virtualenv --python=[path to specific python installation] [path to new virtual environment]
# Activate Virtual Environment:
# source "[virtual environment path and name]/bin/activate"
# Calls to Python and its packages will use the Virtual Environment's
# installations while active.
# Deactivate Virtual Environment:
# deactivate

# Delete virtual environment.
# rm -rf [virtual environment path and name] # remove virtual environment

"${path_python_2718}/bin/python2" -m pip install --user --upgrade virtualenv
"${path_python_2718}/bin/python2" -m virtualenv --version # virtualenv 20.4.7, TCW, 7 July 2021

# Virtual Environment: "ldsc"
# Satisfy dependencies for ldsc package.
# https://github.com/bulik/ldsc/blob/master/environment.yml
"${path_python_2718}/bin/python2" -m virtualenv --python="${path_python_2718}/bin/python2" $path_environment_ldsc

source "${path_environment_ldsc}/bin/activate"
which python2 # "${path_environment_ldsc}/bin/python2" TCW, 7 July 2021
python2 -m pip install bitarray==0.8.3
python2 -m pip install nose==1.3.7
python2 -m pip install pybedtools==0.7.10
python2 -m pip install scipy==0.18.1
python2 -m pip install numpy==1.16.6
python2 -m pip install pandas==0.20.3
python2 -m pip list # TCW, 7 July 2021
deactivate
which python2 # "/usr/bin/python2" TCW, 7 July 2021

################################################################################
# Stuff that still needs to be cleaned up...

path_python="/usr/bin/python3"
path_python_library="/home/tcameronwaller/project_local/python_library"

# Echo each command to console.
set -x

cd ~

# TODO: experiment with python virtual environments...


# Python.
python3 --version
#wget "https://www.python.org/ftp/python/3.9.4/Python-3.9.4.tgz"
#sudo apt update
#sudo apt install python3.9

# TODO: consider using a "requirements.txt" file to specify package names and versions

# Pip.
#python3 -m pip --version
#python3 -m ensurepip --upgrade
#sudo apt install python3-pip
#python3 -m pip3 install --upgrade pip3
pip3 --version

# Install Python packages to specific location.
# pip install --target=[path to installation] --upgrade [package]
# pip install --install-option="--prefix=[path to installation]" --upgrade [package]



# Install specific versions of Python packages
# python -m pip install [requested package]=2.18.4

pip3 install --target=${path_python_library} --upgrade numpy
pip3 install --target=${path_python_library} --upgrade scipy
pip3 install --target=${path_python_library} --upgrade testresources
pip3 install --target=${path_python_library} --upgrade pandas
pip3 install --target=${path_python_library} --upgrade sklearn
pip3 install --target=${path_python_library} --upgrade statsmodels
pip3 install --target=${path_python_library} --upgrade matplotlib
pip3 install --target=${path_python_library} --upgrade networkx

# Uninstall Python packages from specific location.
# A reciprocal uninstall implementation does not yet exist.
# pip uninstall --target=[path] [package]
