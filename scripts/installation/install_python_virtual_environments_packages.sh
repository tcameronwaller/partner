#!/bin/bash
#chmod u+x script.sh

################################################################################
# Paths

# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_python="${path_tools}/python"
path_python_3124="${path_python}/python-3.12.4"
path_python_3111="${path_python}/python_3.11.1"
path_python_396="${path_python}/python_3.9.6"
path_python_3816="${path_python}/python_3.8.16"
path_python_2718="${path_python}/python_2.7.18"

# Python 3 environments.
path_environment_main="${path_tools}/python/environments/main"
path_environment_prs_cs="${path_tools}/python/environments/prs_cs"
path_environment_crossmap="${path_tools}/python/environments/crossmap"
path_environment_gwas2vcf="${path_tools}/python/environments/gwas2vcf"
path_environment_sumstatsrehab="${path_tools}/python/environments/sumstatsrehab"

# Python 2 environments.
path_environment_ldsc="${path_tools}/python/environments/ldsc"

# Initialize directories.
mkdir -p $path_python
mkdir -p $path_python_3124
mkdir -p $path_python_3111
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
# Note: TCW; 23 January 2023
# I do not have root (sudo) privileges on the server and cluster of my current
# organization.
sudo apt update
sudo apt install gcc g++ make wget
sudo apt install dpkg-dev build-essential
sudo apt install zlib1g-dev libbz2-dev libssl-dev uuid-dev libffi-dev libreadline-dev libsqlite3-dev tk-dev libbz2-dev libncurses5-dev libreadline6-dev libgdbm-dev liblzma-dev
sudo apt install libgdbm-compat-dev
sudo apt install libffi-dev libssl-dev openssl-devel
sudo apt install python-dev python3-dev


# Latest
##########
# Python v3.12.4
# date, release: 6 June 2024
# description: https://www.python.org/downloads/release/python-3124/
# installation:
# - system: NCSA, mForge, endocrinology workspace
# - date, installation: 12 July 2024
# - version, installation: v3.12.4
# installation:
# - system: halyard
# - date, installation: 17 July 2024
# - version, installation: v3.12.4
# documentation: https://docs.python.org/3/using/unix.html#getting-and-installing-the-latest-version-of-python
cd "${path_tool}/python"
wget https://www.python.org/ftp/python/3.12.4/Python-3.12.4.tar.xz
tar -xJvf Python-3.12.4.tar.xz # Extract from tar ball with "xz" compression.
rm Python-3.12.4.tar.xz
mv Python-3.12.4 python-3.12.4
cd python-3.12.4
pwd
# Argument "--prefix" needs a full, explicit file path.
# This argument keeps installation within the specific directory path.
./configure --prefix="/.../tool/python/python-3.12.4" # requires full, absolute path to directory
make
make test
# mforge: 468 tests OK.
# mforge: 2 tests failed: test_pathlib test_zipfile
# halyard: 471 tests OK.
# halyard: Result: SUCCESS
make install # mforge: 12 July 2024; halyard: 17 July 2024
"${path_tool}/python/python-3.12.4/bin/python3" -V # "Python 3.12.4"; TCW; 12 July 2024
# Remove installation.
#rm -rf $path_python_3124



##########
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



##########
# Pip
# Install pip for specific installations of Python.
# Note: This execution does not upgrade pip within specific Python installation.
# Note: Hence it is necessary to upgrade within each virtual environment below.
# Note: Activation of the virtual environment changes the "PATH" variable.
# Note: This "PATH" tells Python where to look for the "bin" directory for
# Note: installations.
#"${path_python_396}/bin/python3" -m pip install --user --upgrade pip
#"${path_python_396}/bin/python3" -m pip --version # pip 21.2.1, TCW, 28 July 2021
#"${path_python_2718}/bin/python2" -m pip install --user --upgrade pip
#"${path_python_2718}/bin/python2" -m pip --version # pip 20.3.4, TCW, 6 July 2021



################################################################################
# Python virtual environments
# If you activate a Virtual Environment before installing Python packages to a
# specific path, then the Virtual Environment will remember those paths within
# its own "PYTHONPATH" variable.



##########
# Python3 virtual environment: General Notes
# Python standard module "venv" is the preferred Virtual Environment manager for
# Python version 3.3 or higher.
# https://docs.python.org/3/tutorial/venv.html
# [path to specific python3 installation] -m venv [path to new virtual environment and its name to create]
# source [virtual environment name]/bin/activate
# deactivate
# Delete virtual environment.
# rm -rf [virtual environment path and name] # remove virtual environment



##########
# Python2 virtual environment: General Notes
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


##########
# Python3 environment: "main"
# installation:
# - system: halyard
# - date, installation: 17 July 2024
# - version, installation: v3.12.4
# Initialize virtual environment for Python v3.12.4.
# path to python installation: path_python_3124="${path_tool}/python/python-3.12.4"
# path to environment: path_environment_htseq="${path_tool}/python/environments/main"
# Notice that it is unnecessary to create the last directory ("main") before
# creating the virtual environment.
"${path_python_396}/bin/python3" -m venv --help
"${path_python_396}/bin/python3" -m venv $path_environment_main # TCW; 17 July 2024
# Activate virtual environment.
source "${path_environment_main}/bin/activate"
which python3 # "${path_tool}/python/environments/main/bin/python3"; TCW; 17 July 2024
# Pip installation within virtual environment should not require "sudo" root user permissions.
# Documentation: https://docs.python.org/3/installing/index.html
python3 -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org pip setuptools
python3 -m pip --version # pip 24.0
python3 -m pip install --upgrade pip
python3 -m pip --version # pip 24.1.2
python3 -m pip install --upgrade pytest
python3 -m pip install --upgrade testresources
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade scipy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scikit-learn
python3 -m pip install --upgrade statsmodels
python3 -m pip install --upgrade networkx
python3 -m pip install --upgrade matplotlib
deactivate
which python3


##########
# Python3 environment: "htseq"
# source for HTSeq: https://github.com/htseq/htseq
# documentation for HTSeq: https://htseq.readthedocs.io/en/latest/index.html
# installation:
# - system: NCSA, mForge, endocrinology workspace
# - date, installation: 12 July 2024
# - version, installation: v2.0.8
# Initialize virtual environment for Python v3.12.4.
# path to python installation: path_python_3124="${path_tool}/python/python-3.12.4"
# path to environment: path_environment_htseq="${path_tool}/python/environments/htseq"
# Notice that it is unnecessary to create the last directory ("htseq") before
# creating the virtual environment.
"${path_python_3124}/bin/python3" -m venv --help
"${path_python_3124}/bin/python3" -m venv $path_environment_htseq # TCW; 12 July 2024
source "${path_environment_htseq}/bin/activate"
which python3 # "${path_tool}/python/environments/htseq/bin/python3"; TCW; 12 July 2024
python3 -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org pip setuptools
python3 -m pip install --upgrade pip # TCW; 12 July 2024
python3 -m pip --version # "pip 24.1.2 from ${path_environment_htseq}/lib/python3.12/site-packages/pip (python 3.12)"; TCW; 12 July 2024
python3 -m pip install --upgrade Cython # "Cython-3.0.10"; TCW; 12 July 2024
python3 -m pip install --upgrade pytest # "dependencies ... pytest-8.2.2"; TCW; 12 July 2024
python3 -m pip install --upgrade numpy # "numpy-2.0.0"; TCW; 12 July 2024
python3 -m pip install --upgrade scipy # "scipy-1.14.0"; TCW; 12 July 2024
python3 -m pip install --upgrade pandas # "dependencies ... pandas-2.2.2"; TCW; 12 July 2024
python3 -m pip install --upgrade matplotlib # "dependencies ... matplotlib-3.9.1"; TCW; 12 July 2024
python3 -m pip install --upgrade pysam # "pysam-0.22.1"; TCW; 12 July 2024
python3 -m pip install --upgrade anndata # "dependencies ... anndata-0.10.8"; TCW; 12 July 2024
python3 -m pip install --upgrade loompy # "dependencies ... loompy-3.0.7"; TCW; 12 July 2024
python3 -m pip install --upgrade pyBigWig # "pyBigWig-0.3.23"; TCW; 12 July 2024
python3 -m pip install --upgrade swig # "swig-4.2.1"; TCW; 12 July 2024
python3 -m pip install --upgrade HTSeq # "2.0.8"; TCW; 12 July 2024
python3 -m HTSeq.scripts.count --version # "2.0.8"; TCW; 12 July 2024
python3 -m HTSeq.scripts.count --help
deactivate
which python3



##########
# Python3 environment: "crossmap"
# Documentation for Python 2 CrossMap: "https://pythonhosted.org/CrossMap/"
# Documentation for Python 3 CrossMap: "https://sourceforge.net/projects/crossmap/"
# Documentation for Python 3 CrossMap: "http://crossmap.sourceforge.net/"
# Initialize Python 3.9.6 virtual environment.
# Installation of Python 3.9.6: path_python_396="${path_python}/python_3.9.6"
# Path to environment: path_environment_crossmap="${path_tools}/python/environments/crossmap"
"${path_python_396}/bin/python3" -m venv $path_environment_crossmap # TCW; 10 March 2023
source "${path_environment_crossmap}/bin/activate"
which python3 # "${path_environment_crossmap}/bin/python3"; TCW; 10 March 2023
python3 -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org pip setuptools
python3 -m pip install --upgrade pip # TCW; 10 March 2023
# Successfully installed pip-23.0.1
python3 -m pip --version # "pip 23.0.1 from ${path_environment_crossmap}/lib/python3.9/site-packages/pip (python 3.9)"; TCW; 10 March 2023
python3 -m pip install --upgrade CrossMap # TCW; 10 March 2023
# "Successfully installed CrossMap-0.6.5 bx-python-0.9.0 cython-0.29.33 numpy-1.24.2 pyBigWig-0.3.20 pysam-0.20.0"
CrossMap.py --help # call "CrossMap.py" directly as an executable
CrossMap.py bed --help
CrossMap.py vcf --help
deactivate
which python3



##########
# Python3 environment: "gwas2vcf"
# Host of GWAS2VCF: https://github.com/MRCIEU/gwas2vcf
# Documentation for GWAS2VCF: https://mrcieu.github.io/gwas2vcf/install/#dbsnp
# Package GWAS2VCF requires Python 3.8.
# Initialize Python 3.8 virtual environment.
# Installation of Python 3.8: path_python_3816="${path_python}/python_3.8.16"
# Path to environment: path_environment_gwas2vcf="${path_tools}/python/environments/gwas2vcf"
"${path_python_3816}/bin/python3" -m venv $path_environment_gwas2vcf # TCW; 23 January 2023
source "${path_environment_gwas2vcf}/bin/activate"
which python3 # "${path_environment_gwas2vcf}/bin/python3"; TCW; 23 January 2023
python3 -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org pip setuptools
python3 -m pip --version # "pip 22.0.4 from ${path_environment_gwas2vcf}/lib/python3.8/site-packages/pip (python 3.8)"; TCW; 23 January 2023
python3 -m pip install --upgrade pip
python3 -m pip --version # "pip 22.3.1 from ${path_environment_gwas2vcf}/lib/python3.8/site-packages/pip (python 3.8)"; TCW; 23 January 2023
# Dependencies or requirements of GWAS2VCF package.
# https://github.com/MRCIEU/gwas2vcf/blob/master/requirements.txt
python3 -m pip install pytest==5.3.5 # TCW; 23 January 2023
python3 -m pip install biopython==1.72 # TCW; 23 January 2023
python3 -m pip install pysam==0.15.2 # TCW; 23 January 2023
python3 -m pip install marshmallow==2.18.1 # TCW; 23 January 2023
python3 -m pip install numpy==1.22.0 # TCW; 23 January 2023
python3 -m pip install Cython==0.29.17 # TCW; 23 January 2023
# Example: pip install https://git+github.com/<owner_name>/<repo_name>.git.@<version#>
python3 -m pip install git+https://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph # TCW; 23 January 2023
deactivate
which python3



##########
# Python3 environment: "sumstatsrehab"



################################################################################


##########

# repository and documentation: https://github.com/getian107/PRScs/tree/master

# Python3 virtual environment: "prs_cs"
# path_python_396: full file path to installation of Python 3.9.6 or alternative.
# path_environment_prs_cs: full file path for installation of virtual environment.

# Probably best to specify versions of packages to install using syntax below.
#python3 -m pip install numpy==1.22.3 # Or more recent version that works.

# Python3 virtual environment: "prs_cs"
"${path_python_396}/bin/python3" -m venv --help # TCW; 28 April 2022
# ./python/python_3.9.6/bin/python3 -m venv ./python/environments/prs_cs # TCW; 28 April 2022
"${path_python_396}/bin/python3" -m venv $path_environment_prs_cs # TCW; 28 April 2022
# source ./python/environments/prs_cs/bin/activate # TCW; 28 April 2022
source "${path_environment_prs_cs}/bin/activate"
which python3 # "${path_environment_main}/bin/python3" TCW, 28 July 2021
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



##########
# Python2 virtual environment: "ldsc"
"${path_python_2718}/bin/python2" -m pip install --user --upgrade pip
"${path_python_2718}/bin/python2" -m pip --version # pip 20.3.4, TCW, 6 July 2021
"${path_python_2718}/bin/python2" -m pip install --user --upgrade virtualenv
"${path_python_2718}/bin/python2" -m virtualenv --version # virtualenv 20.4.7, TCW, 7 July 2021
# Virtual Environment: "ldsc"
# Satisfy dependencies for ldsc package.
# https://github.com/bulik/ldsc/blob/master/requirements.txt
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
