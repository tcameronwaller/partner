#!/bin/bash

#chmod u+x script.sh

# Echo each command to console.
set -x

# Set working directory.
cd ~



##################################################
# Install tools for local personal and professional use
##################################################

# Install general tools.
# "software-properties-common" increases "apt" functionality for commands such as "add-apt-repository".
sudo apt install software-properties-common
sudo apt install wget
sudo apt install curl
sudo apt install unzip

# Install Tweak tool to customize Ubuntu's appearance and behavior.
sudo add-apt-repository universe
sudo apt install gnome-tweaks


##########
# Libre Office
# site: https://www.libreoffice.org/
# note:
# While Libre Office is available as a Snap application on Ubuntu, this version
# tends to have performance issues of lag.
# uninstall snap version
# sudo snap remove --purge libreoffice
# installation:
# - date, installation: 5 March 2025
# - version, installation: v25.2.1
# confirm installation of dependency, Java Runtime Environment
java -version
wget https://www.libreoffice.org/donate/dl/deb-x86_64/25.2.1/en-US/LibreOffice_25.2.1_Linux_x86-64_deb.tar.gz
wget https://www.libreoffice.org/donate/dl/deb-x86_64/25.2.1/en-US/LibreOffice_25.2.1_Linux_x86-64_deb.tar.gz
tar -xzvf ./LibreOffice_25.2.1_Linux_x86-64_deb.tar.gz -C ./
#tar --extract --file="LibreOffice_25.2.1_Linux_x86-64_deb.tar.gz" --gzip --verbose
cd ./LibreOffice_25.2.1_Linux_x86-64_deb/DEBS/
sudo dpkg -i ./LibreOffice_25.2.1_Linux_x86-64_deb/DEBS/*.deb
# uninstall
sudo dpkg -r packagename


##########
# Atom editor for text and code
# site: https://atom-editor.cc/
# source: https://github.com/atom/atom/releases/tag/v1.60.0
# version: v1.60.0
# date, release: 3 March 2023
# date, installation: 10 December 2024
# note:
# - The Atom Editor project was formally discontinued in 2022 or 2023.
# - Install from archive.
# download
$ wget https://github.com/atom/atom/releases/download/v1.60.0/atom-amd64.deb
# install
$ sudo apt install ./atom-amd64.deb
# or alternatively
$ sudo dpkg -i ./atom-amd64.deb
# uninstall
$ sudo apt remove atom
$ sudo apt autoremove
# Notes from previous installation while repository was active.
#curl -sL https://packagecloud.io/AtomEditor/atom/gpgkey | sudo apt-key add -
#sudo sh -c 'echo "deb [arch=amd64] https://packagecloud.io/AtomEditor/atom/any/ any main" > /etc/apt/sources.list.d/atom.list'
#wget -qO - https://packagecloud.io/AtomEditor/atom/gpgkey | sudo apt-key add -
#sudo add-apt-repository "deb [arch=amd64] https://packagecloud.io/AtomEditor/atom/any/ any main"
#sudo apt update
#sudo apt install atom


# Discord
# Or use the web browser application.
cd ~/Downloads
wget https://dl.discordapp.net/apps/linux/0.0.40/discord-0.0.49.deb
#wget https://discord.com/api/download?platform=linux&format=deb
#wget "https://discord.com/api/download?platform=linux&format=deb"
sudo apt install ./discord-0.0.63.deb # TCW; 9 August 2024
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
wget https://zoom.us/client/6.1.0.198/zoom_amd64.deb
#wget https://zoom.us/client/5.15.11.7239/zoom_amd64.deb
#wget https://zoom.us/client/5.17.10.3512/zoom_amd64.deb
#wget https://zoom.us/client/latest/zoom_amd64.deb
#wget https://zoom.us/linux/download/pubkey
#rpm --import package-signing-key.pub
#rpm --import pubkey
sudo dpkg -i ./zoom_amd64.deb
#sudo apt install ./zoom_amd64.deb
sudo apt update
sudo apt upgrade
sudo apt remove zoom

# darktable
# description: edit digital photographs in raw format
# site: https://www.darktable.org/
# Install the "universal package formats" using flatpak.
flatpak install flathub darktable

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
# review: TCW; 5 February 2025
sudo add-apt-repository ppa:inkscape.dev/stable
sudo apt update
sudo apt install inkscape
sudo apt upgrade
sudo apt remove inkscape

# Audacity
# site: https://www.audacityteam.org/
# option: install from Ubuntu App Store
sudo snap connect audacity:alsa
# "If you choose to continue and later wish to change your mind, you can remove the file
# /home/tcameronwaller/snap/audacity/common/Always use PulseAudio."
#cd ~/Downloads
#wget https://github.com/audacity/audacity/releases/download/Audacity-3.0.3/audacity-linux-3.0.3-x86_64.AppImage
#chmod +x ./audacity-linux-3.0.3-x86_64.AppImage
#./audacity-linux-3.0.3-x86_64.AppImage

# VLC Media Player
# site: https://www.videolan.org/vlc/download-ubuntu.html
# date, installation: 16 February 2025
sudo snap install vlc
sudo snap remove --purge vlc

# Rhythmbox
# site: http://www.rhythmbox.org/


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
# date: 15 February 2025
cd /folder/where/you/saved/the/file
# Make file executable as a program.
chmod +x ./Ardour-8.11.0-x86_64.run
#/bin/sh ./Ardour-8.11.0-x86_64.run
./Ardour-8.11.0-x86_64.run
# !!! WARNING !!! - Your system seems to use frequency scaling.
# This can have a serious impact on audio latency.
# For best results turn it off, e.g. by choosing the 'performance' governor.
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



#
