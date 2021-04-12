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

# Discord
cd ~/Downloads
wget https://dl.discordapp.net/apps/linux/0.0.14/discord-0.0.14.deb
#wget "https://discord.com/api/download?platform=linux&format=deb"
sudo apt install ./discord-0.0.14.deb
sudo apt remove discord

# Zoom
# https://zoom.us/download
# - Ubuntu Linux
# - 64 bit
# - 16.04+
cd ~/Downloads
wget https://zoom.us/client/latest/zoom_amd64.deb
sudo apt install ./zoom_amd64.deb
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
sudo add-apt-repository ppa:inkscape.dev/stable
sudo apt update
sudo apt install inkscape
