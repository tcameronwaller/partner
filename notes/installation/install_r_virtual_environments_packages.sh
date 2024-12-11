

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tools=$(<"$path_directory_paths/path_directory_tools.txt")
path_directory_repository_partner=$(<"$path_directory_paths/path_directory_repository_partner.txt")
path_directory_scripts="$path_directory_repository_partner/scripts"


###############################################################################
# Installation of R
# site: https://www.r-project.org/
# CRAN mirror host: Iowa State University in Ames, Iowa
# site: https://mirror.las.iastate.edu/CRAN/

##########
# R v4.4.1
# description: https://mirror.las.iastate.edu/CRAN/doc/manuals/r-release/NEWS.html
# version: v4.4.2
# date, release: 31 October 2024
# installation:
# - installation from source
# - system: halyard
# - date, installation: 10 December 2024
# - version, installation: v4.4.2

# Install dependencies.
# https://docs.posit.co/resources/install-r-source.html
# https://cran.r-project.org/bin/linux/ubuntu/
#sudo sed -i.bak "/^#.*deb-src.*universe$/s/^# //g" /etc/apt/sources.list # <-- do this in similar way to how I did it for Python
sudo apt install --no-install-recommends software-properties-common dirmngr

# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo add-apt-repository --enable-source https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/"

sudo nano /etc/apt/sources.list
# "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/"
# "deb-src https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/"
less /etc/apt/sources.list

sudo apt update
sudo apt install r-base-dev
sudo apt build-dep r-base
sudo apt update
# Install R itself.
cd "${path_directory_tools}/r"
wget https://mirror.las.iastate.edu/CRAN/src/base/R-4/R-4.4.2.tar.gz
tar -xzvf R-4.4.2.tar.gz # Extract from tar ball with GZip compression.
rm R-4.4.2.tar.gz
mv R-4.4.2 r-4.4.2
cd r-4.4.2
pwd
# Argument "--prefix" needs a full, explicit file path.
# This argument keeps installation within the specific directory path.
./configure \
--prefix="/.../tool/r/r-4.4.2" \
--enable-R-shlib \
--enable-memory-profiling
make
make check # synonym for "make test"?; halyard; 10 December 2024
make install # halyard; 10 December 2024
# Rscript.
"${path_directory_tools}/r/r-4.4.1/bin/Rscript" --help
"${path_directory_tools}/r/r-4.4.1/bin/Rscript" --version # "Rscript (R) version 4.4.2 (2024-10-31)"; halyard; 10 December 2024
# Initialize R console.
# Execution from Bash terminal.
# user@machine:/path$ "${path_tool}/r/r-4.4.1/bin/R"
# Execution from R console.
# > help.start()
# > ?mean
# Print path to library in which R installs packages.
# > getOption("repos") # Print repositories from which R will install packages.
# > .libPaths() # "$path_tool/r/r-4.4.1/lib/R/library"
# > library()
# > library(lib.loc='~/r-packages')
# > lapply(.libPaths(), list.files) # Print packages available in each environmental library.
# > q()
# Remove installation.
#rm -rf $path_r_441



###############################################################################
# Installation of Bioconductor repository and other packages

# Initialize R console and execute installations interactively.
# TCW; 11 December 2024

# Execution from Bash terminal.
"${path_directory_tools}/r/r-4.4.1/bin/Rscript" "$path_directory_scripts/r/install_packages.R"

#Rscript -e 'renv::run("/path/to/myscript.R")'


# End
