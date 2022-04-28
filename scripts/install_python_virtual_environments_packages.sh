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

path_environment_prs_cs="${path_tools}/python/environments/prs_cs"

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

# Install dependencies before installing Python2 and Python3 from source.
# https://tttapa.github.io/Pages/Ubuntu/Software-Installation/Python.html
sudo apt install dpkg-dev build-essential
sudo apt install zlib1g-dev libbz2-dev libssl-dev uuid-dev libffi-dev libreadline-dev libsqlite3-dev tk-dev libbz2-dev libncurses5-dev libreadline6-dev libgdbm-dev liblzma-dev
sudo apt install libgdbm-compat-dev
sudo apt install libffi-dev libssl-dev openssl-devel
sudo apt install python-dev python3-dev

# Latest
# Python 3.9.6, release date: 28 June 2021
# https://www.python.org/downloads/release/python-396/
# https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz
cd $path_python_396
wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz # TCW, 6 July 2021
tar -xzvf Python-3.9.6.tgz
#find ~/python -type d | xargs chmod 0755
cd Python-3.9.6
# Argument "--prefix" needs a full, explicit path.
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

# Note: this doesn't upgrade pip in the specific python3 installation
# Note: hence necessary to upgrade within the environment
"${path_python_396}/bin/python3" -m pip install --user --upgrade pip
"${path_python_396}/bin/python3" -m pip --version # pip 21.2.1, TCW, 28 July 2021

"${path_python_2718}/bin/python2" -m pip install --user --upgrade pip
"${path_python_2718}/bin/python2" -m pip --version # pip 20.3.4, TCW, 6 July 2021

################################################################################
# Python virtual environments
# If you activate a Virtual Environment before installing Python packages to a
# specific path, then the Virtual Environment will remember those paths within
# its own "PYTHONPATH" variable.

################################################################################
# Python3 virtual environment
# Python standard module "venv" is the preferred Virtual Environment manager for
# Python version 3.3 or higher.
# https://docs.python.org/3/tutorial/venv.html
# [path to specific python3 installation] -m venv [path to new virtual environment and its name to create]
# source [virtual environment name]/bin/activate
# deactivate

"${path_python_396}/bin/python3" -m venv --help # TCW, 28 July 2021

"${path_python_396}/bin/python3" -m venv $path_environment_main # TCW, 28 July 2021

# Activate virtual environment.
# Unable to install with "--user" flag within virtual environment.
# ERROR: Can not perform a '--user' install. User site-packages are not visible in this virtualenv.
source "${path_environment_main}/bin/activate"
which python3 # "${path_environment_main}/bin/python3" TCW, 28 July 2021
# Pip installation within virtual environment should not require "sudo" permissions.
python3 -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org pip setuptools
python3 -m pip --version # "pip 21.1.3 from ${path_environment_main}/lib/python3.9/site-packages/pip (python 3.9)" TCW, 28 July 2021
python3 -m pip install --upgrade pip
python3 -m pip --version # "pip 22.0.4 from ${path_environment_main}/lib/python3.9/site-packages/pip (python 3.9)" TCW, 12 April 2022
python3 -m pip install --upgrade numpy # "numpy-1.21.1" TCW, 28 July 2021
python3 -m pip install --upgrade scipy # "scipy-1.7.0" TCW, 28 July 2021
python3 -m pip install --upgrade testresources # "pbr-5.6.0 testresources-2.0.1" TCW, 28 July 2021
python3 -m pip install --upgrade pandas # "pandas-1.3.1 python-dateutil-2.8.2 pytz-2021.1 six-1.16.0" TCW, 28 July 2021
python3 -m pip install --upgrade sklearn # "joblib-1.0.1 scikit-learn-0.24.2 sklearn-0.0 threadpoolctl-2.2.0" TCW, 28 July 2021
python3 -m pip install --upgrade statsmodels # "patsy-0.5.1 statsmodels-0.12.2" TCW, 28 July 2021
python3 -m pip install --upgrade matplotlib # "cycler-0.10.0 kiwisolver-1.3.1 matplotlib-3.5.1 pillow-8.3.1 pyparsing-2.4.7" TCW, 12 April 2022
python3 -m pip install --upgrade networkx # "networkx-2.6.2" TCW, 28 July 2021
deactivate
which python3

# Delete virtual environment.
# rm -rf [virtual environment path and name] # remove virtual environment

################################################################################
# Python3 virtual environment
# Python standard module "venv" is the preferred Virtual Environment manager for
# Python version 3.3 or higher.
# https://docs.python.org/3/tutorial/venv.html
# [path to specific python3 installation] -m venv [path to new virtual environment and its name to create]
# source [virtual environment name]/bin/activate
# deactivate

# ./python/python_3.9.6/bin/python3 -m venv --help
"${path_python_396}/bin/python3" -m venv --help # TCW; 28 April 2022

# ./python/python_3.9.6/bin/python3 -m venv ./python/environments/prs_cs # TCW; 28 April 2022
"${path_python_396}/bin/python3" -m venv $path_environment_prs_cs # TCW; 28 April 2022

# Activate virtual environment.
# Unable to install with "--user" flag within virtual environment.
# ERROR: Can not perform a '--user' install. User site-packages are not visible in this virtualenv.

# source ./python/environments/prs_cs/bin/activate # TCW; 28 April 2022
source "${path_environment_prs_cs}/bin/activate"
which python3 # "${path_environment_main}/bin/python3" TCW, 28 July 2021
# Pip installation within virtual environment should not require "sudo" permissions.
python3 -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org pip setuptools
python3 -m pip --version # "pip 21.1.3 from ${path_environment_main}/lib/python3.9/site-packages/pip (python 3.9)" TCW; 28 April 2022
python3 -m pip install --upgrade pip # TCW; 28 April 2022
python3 -m pip --version # "pip 22.0.4 from ${path_environment_main}/lib/python3.9/site-packages/pip (python 3.9)" TCW; 28 April 2022
python3 -m pip install --upgrade numpy # "numpy-1.22.3" TCW; 28 April 2022
python3 -m pip install --upgrade scipy # "scipy-1.8.0" TCW; 28 April 2022
python3 -m pip install --upgrade h5py # "h5py-3.6.0" TCW; 28 April 2022
python3 -m pip install --upgrade testresources # "pbr-5.8.1 testresources-2.0.1" TCW; 28 April 2022
deactivate
which python3

# Delete virtual environment.
# rm -rf [virtual environment path and name] # remove virtual environment


################################################################################
# Python2 virtual environment

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
################################################################################
################################################################################
# Notes

# Install Python packages to specific location.
# Prioritize packages at this path by appending it to the "PYTHONPATH" global variable.
#path_python_library="/.../python_library"
#PYTHONPATH=${path_python_library}:$PYTHONPATH
#export PYTHONPATH
# pip install --target=[path to installation] --upgrade [package]
# pip install --install-option="--prefix=[path to installation]" --upgrade [package]

# Install specific versions of Python packages
# python -m pip install [requested package]=2.18.4

# Uninstall Python packages from specific location.
# A reciprocal uninstall implementation does not yet exist.
# pip uninstall --target=[path] [package]
