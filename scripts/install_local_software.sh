#!/bin/bash

#chmod u+x script.sh

# Echo each command to console.
set -x

# Set working directory.
cd ~

# Install general tools.
# "software-properties-common" increases "apt" functionality for commands such as "add-apt-repository".
sudo apt install software-properties-common
sudo apt install wget
sudo apt install curl

# Install Tweak tool to customize Ubuntu's appearance and behavior.
sudo add-apt-repository universe
sudo apt install gnome-tweaks

# Install Atom syntax code editor.
#curl -sL https://packagecloud.io/AtomEditor/atom/gpgkey | sudo apt-key add -
#sudo sh -c 'echo "deb [arch=amd64] https://packagecloud.io/AtomEditor/atom/any/ any main" > /etc/apt/sources.list.d/atom.list'
wget -qO - https://packagecloud.io/AtomEditor/atom/gpgkey | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://packagecloud.io/AtomEditor/atom/any/ any main"
sudo apt update
sudo apt install atom

# Cytoscape
# site: https://cytoscape.org/index.html
# last installation: TCW; 22 February 2024
# 1. Install Java 17
apt update
apt upgrade
apt install openjdk-17-jdk openjdk-17-jre
java -version
# 2. Install Cytoscape
cd ~/Downloads
wget https://github.com/cytoscape/cytoscape/releases/download/3.10.1/Cytoscape_3_10_1_unix.sh
bash ~/Downloads/Cytoscape_3_10_1_unix.sh
# installation path: /home/tcameronwaller/Cytoscape_v3.10.1

# Discord
# Or use the web browser application.
cd ~/Downloads
wget https://dl.discordapp.net/apps/linux/0.0.40/discord-0.0.49.deb
#wget https://discord.com/api/download?platform=linux&format=deb
#wget "https://discord.com/api/download?platform=linux&format=deb"
sudo apt install ./discord-0.0.49.deb
sudo apt remove discord

# Spotify
# https://www.spotify.com/de-en/download/linux/
# version: spotify 1.2.31.1205.g4d59ad7c
# date, installation: 26 March 2024
snap install spotify

# Zoom
# "W: Repository is broken: zoom:amd64 (= 5.7.31792.0820) has no Size information"
# https://zoom.us/download
# - Ubuntu Linux
# - 64 bit
# - 16.04+
sudo apt install libglib2.0-0 libxcb-shape0 libxcb-shm0 libxcb-xfixes0 libxcb-randr0 libxcb-image0 libfontconfig1 libgl1-mesa-glx libxi6 libsm6 libxrender1 libpulse0 libxcomposite1 libxslt1.1 libsqlite3-0 libxcb-keysyms1 libxcb-xtest0 ibus
sudo apt install libxcb-cursor0 # TCW; 23 March 2023
# libgstreamer-plugins-base0.10-0
cd ~/Downloads
#wget https://zoom.us/client/5.15.11.7239/zoom_amd64.deb
wget https://zoom.us/client/5.17.10.3512/zoom_amd64.deb
#wget https://zoom.us/client/latest/zoom_amd64.deb
#wget https://zoom.us/linux/download/pubkey
#rpm --import package-signing-key.pub
#rpm --import pubkey
sudo dpkg -i ./zoom_amd64.deb
#sudo apt install ./zoom_amd64.deb
sudo apt update
sudo apt upgrade
sudo apt remove zoom

# GIMP
#sudo add-apt-repository ppa:otto-kesselgulasch/gimp
#sudo apt update
#sudo apt install gimp
sudo apt install flatpak
flatpak install https://flathub.org/repo/appstream/org.gimp.GIMP.flatpakref
flatpak update
flatpak run org.gimp.GIMP//stable

# PulseEffects
# Audio effects
sudo add-apt-repository ppa:mikhailnov/pulseeffects
sudo apt update
sudo apt install pulseaudio pulseeffects --install-recommends
# sudo apt remove pulseeffects
# or use FlatPak
flatpak install flathub com.github.wwmm.pulseeffects
flatpak run com.github.wwmm.pulseeffects

# Tweaks
# Customize Ubuntu user interface themes and additional system controls.
# Change suspend behavior on closing computer lid.
sudo add-apt-repository universe
sudo apt install gnome-tweak-tool

# Inkscape
# website: https://inkscape.org/
# review: TCW; 25 October 2022
sudo add-apt-repository ppa:inkscape.dev/stable
sudo apt update
sudo apt install inkscape
sudo apt remove inkscape

# Audacity
# https://www.audacityteam.org/
cd ~/Downloads
wget https://github.com/audacity/audacity/releases/download/Audacity-3.0.3/audacity-linux-3.0.3-x86_64.AppImage
chmod +x ./audacity-linux-3.0.3-x86_64.AppImage
./audacity-linux-3.0.3-x86_64.AppImage

# MuseScore
# https://musescore.org/en
#cd ~/Downloads
#wget https://musescore.org/en/download/musescore-x86_64.AppImage
#chmod +x
#./musescore-x86_64.AppImage
sudo add-apt-repository ppa:mscore-ubuntu/mscore-stable
sudo apt-get update
sudo apt install musescore

# PulseAudio and PulseAudio Volume Control
# PulseAudio improves performance of Audacity for audio recording.
# PulseAudio might conflict with Jack Audio
sudo apt install pulseaudio
sudo apt install pavucontrol

# Jack Audio Connection Kit
# Utility for high-resolution audio recording without interference from
# central processing unit.
sudo apt install jackd
sudo apt install qjackctl
whereis jackd
whereis qjackctl
qjackctl &
# Enable Realtime scheduling and memory locking explicitly for user.
# Make sure file exists: "/etc/security/limits.d/audio.conf".
# If file exists: "/etc/security/limits.d/audio.conf.disabled" then rename to "/etc/security/limits.d/audio.conf".
# Run command: "sudo usermod -a -G audio tcameronwaller"
# Restart system to enact changes.

# Manager for CPU Frequency Scaling
# Select "performance" setting to optimize CPU performance frequency.
# Important for recording audio with low latency.
sudo apt install linux-tools-common
sudo apt install linux-tools-generic
sudo apt install linux-tools-5.11.0-27-generic
sudo cpupower -c all frequency-set -g performance

# Ardour
# http://ardour.org/first_time_linux.html
cd /folder/where/you/saved/the/file
/bin/sh ./<DOWNLOADED_FILENAME>.run
# Uninstall
cd /opt
/bin/sh ./Ardour_<VERSION>.uninstall.sh

# Linux Multi-Media Studio (LMMS)
# https://lmms.io
# Set the download file to executable.
# Run the AppImage directly.
cd ~/Downloads
wget https://github.com/LMMS/lmms/releases/download/v1.2.2/lmms-1.2.2-linux-x86_64.AppImage
chmod +x ./lmms-1.2.2-linux-x86_64.AppImage
./lmms-1.2.2-linux-x86_64.AppImage

##########
# CHIRP
# review: TCW; 15 September 2022
# https://chirp.danplanet.com/projects/chirp
# https://chirp.danplanet.com/projects/chirp/wiki/Running_Under_Linux
# https://chirp.danplanet.com/projects/chirp/wiki/Beginners_Guide
# https://chirp.danplanet.com/projects/chirp/wiki/MemoryEditorColumns
# Installation on Linux Ubuntu via Flatpak repository.
sudo apt install flatpak
flatpak update
sudo flatpak remote-add --if-not-exists flathub https://flathub.org/repo/flathub.flatpakrepo
flatpak update
cd ~/Downloads
wget https://trac.chirp.danplanet.com/chirp_daily/LATEST/chirp-daily-20220911.flatpak
sudo flatpak install chirp-daily-*.flatpak
flatpak run com.danplanet.chirp
# Allow user to access USB ports.
sudo addgroup tcameronwaller dialout
sudo usermod -aG dialout tcameronwaller
dmesg | grep tty # determine to which USB port the radio is connected
# /dev/ttyUSB0

flatpak uninstall com.danplanet.chirp

##########
# FLDigi
# review: TCW; 17 September 2022
# http://www.w1hkj.com/
# https://sourceforge.net/p/fldigi/wiki/debian_howto/
# http://www.w1hkj.com/LaunchpadInstall.html
# https://launchpad.net/~ubuntu-hams-updates/+archive/ubuntu/ppa
# Installation on Linux Ubuntu via PPA Repository.
sudo add-apt-repository ppa:kamalmostafa/fldigi
sudo add-apt-repository ppa:ubuntu-hams-updates/ppa
sudo apt update
sudo apt install fldigi

##########
# AndFlmsg
# review: TCW; 17 September 2022
# This is an application for Android with some functionality similar to FLDigi.
# Follow the instructions for installation on Android.
# It is necessary to download and open the ".apk" file and allow installation
# from a non standard source.
# http://www.w1hkj.com/
# http://www.w1hkj.com/vk2eta/
# http://www.w1hkj.com/files/AndFlmsg/INSTALL.txt
# - - Installation instructions.
# https://sourceforge.net/projects/fldigi/files/AndFlmsg/
# http://www.w1hkj.com/files/AndFlmsg/
# http://www.w1hkj.com/files/AndFlmsg/AndFlmsg_V1.5.0-20210812.apk


##########
# Git, GitHub
# site: https://github.com/
# After renaming a repository, redirect git commits to new repository on GitHub.
git remote set-url origin https://github.com/tcameronwaller/partner.git
git remote set-url origin https://github.com/tcameronwaller/psychiatry_biomarkers.git

##########
# BGZip
# Last installation: TCW; 15 February 2023 (on NCSA server)
# A bioinformatics tool for compression in GZip format with Tabix index.
# http://www.htslib.org/doc/bgzip.html
# The HTSLib includes BGZip and Tabix.
# http://www.htslib.org/download/
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjvf ./htslib-1.16.tar.bz2 -C ./
cd ./htslib-1.16
./configure --prefix=../htslib/1.16 # requires absolute directory path
make
make install
cd ../htslib/1.16/bin/
./bgzip --help



##########
# BCFTools
# Last installation: TCW; 15 February 2023 (on NCSA server)
# A bioinformatic tool for tasks on genotype files in Variant Call Format (VCF).
# Documentation from SamTools (https://samtools.github.io/bcftools/howtos/index.html).
# Documentation manual (https://samtools.github.io/bcftools/bcftools.html).
# https://samtools.github.io/bcftools/
# http://www.htslib.org/download/
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xjvf ./bcftools-1.16.tar.bz2 -C ./
cd ./bcftools-1.16
./configure --prefix=../bcftools/1.16 # requires absolute directory path
make
make install
cd ../bcftools/1.16/bin/
./bcftools --help



##########
# GCTB
# Installation date: TCW; 13 January 2023
# A package of tools for Genome-wide Complex Trait Bayesian (GCTB) analysis
# Includes the tool SBayesR for calculation of Polygenic Scores (PGS).
# Authors provided a "statically linked 64-bit Linux executable, gctb".
# Installation instructions for custom, local compilation (https://cnsgenomics.com/software/gctb/download/README.html).
cd ./ # navigate to directory in which to install program.
wget "https://cnsgenomics.com/software/gctb/download/gctb_2.04.3_Linux.zip"
unzip "./gctb_2.04.3_Linux.zip"
./gctb_2.04.3_Linux/gctb --help

##########
# LDSC
# A package of tools to estimate autosome-wide (excluding sex chromosomes) SNP
# heritability and genetic correlation by Linkage Disequilibrium (LD) Score
# Regression.
# PubMed: ___
# documentation:
# LDSC is a Python package, and it is necessary to run this program within a
# special Python environment.
# Refer to the script "install_python_virtual_environments_packages.sh".
# Navigate to the directory in which to install program.
cd ./ldsc/
# Copy GitHub repository.
git clone https://github.com/bulik/ldsc.git

##########
# PRS-CSX
# Install PRS-CS and PRS-CSX.
cd ./prs_cs/ # Navigate to the directory in which to install program.
git clone https://github.com/getian107/PRScs.git
git clone https://github.com/getian107/PRScsx.git




##########
# GWAS2VCF
# Installation date: TCW; 23 January 2023
# A package of tools to translate GWAS summary statistics to the GWAS-VCF
# format.
# PubMed: 33441155
# GWAS-VCF format specification: https://github.com/MRCIEU/gwas-vcf-specification
# Host of GWAS2VCF: https://github.com/MRCIEU/gwas2vcf
# Documentation for GWAS2VCF: https://mrcieu.github.io/gwas2vcf/install/#dbsnp
# GWAS2VCF is a Python package, and it is necessary to run this program within a
# special Python environment.
# Refer to the script "install_python_virtual_environments_packages.sh".
# Access the GWAS2VCF package.
cd ./tools
mkdir -p ./gwas2vcf
cd ./gwas2vcf # Navigate to the directory in which to install program.
git clone https://github.com/MRCIEU/gwas2vcf.git
python3 ./gwas2vcf/main.py # Execute within a Python virtual environment with dependencies.



##########
# GZ-Sort
# Installation date: TCW; 7 February 2023
# site: http://kmkeen.com/gz-sort/
# repository: https://github.com/keenerd/gz-sort
git clone https://github.com/keenerd/gz-sort; cd gz-sort; make; ./gz-sort -h



##########
# CrossMap
# CrossMap is a Python package, and it is necessary to run this program within a
# special Python environment.
# Refer to the script "install_python_virtual_environments_packages.sh".


##########
# R for Statistical Computing
# https://cran.r-project.org/
# Install R from the CRAN Repository (not the Ubuntu repository).
sudo apt update -qq
sudo apt install --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt install --no-install-recommends r-base
sudo add-apt-repository ppa:c2d4u.team/c2d4u4.0+



#
