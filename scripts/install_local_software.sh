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
# Or use the web browser application.
cd ~/Downloads
wget https://dl.discordapp.net/apps/linux/0.0.15/discord-0.0.15.deb
#wget "https://discord.com/api/download?platform=linux&format=deb"
sudo apt install ./discord-0.0.14.deb
sudo apt remove discord

# Zoom
# "W: Repository is broken: zoom:amd64 (= 5.7.31792.0820) has no Size information"
# https://zoom.us/download
# - Ubuntu Linux
# - 64 bit
# - 16.04+
sudo apt install libglib2.0-0 libxcb-shape0 libxcb-shm0 libxcb-xfixes0 libxcb-randr0 libxcb-image0 libfontconfig1 libgl1-mesa-glx libxi6 libsm6 libxrender1 libpulse0 libxcomposite1 libxslt1.1 libsqlite3-0 libxcb-keysyms1 libxcb-xtest0 ibus
# libgstreamer-plugins-base0.10-0
cd ~/Downloads
wget https://zoom.us/client/latest/zoom_amd64.deb
wget https://zoom.us/linux/download/pubkey
#rpm --import package-signing-key.pub
rpm --import pubkey
sudo dpkg -i ./zoom_amd64.deb
#sudo apt install ./zoom_amd64.deb
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
