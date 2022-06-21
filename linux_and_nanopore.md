
# Install Ubuntu on a Dell desktop configured for the Unified Extensible firmware interface (UEFI) BIOS
https://www.dell.com/support/kbdoc/en-ca/000131655/how-to-install-ubuntu-linux-on-your-dell-pc
## Create a bootable USB flash drive
1. Download Ubuntu desktop image, the version we are using is [Ubuntu 22.04 LTS](https://ubuntu.com/download/desktop)
2. Follow the Canonical ubuntu totorials on how to [Create a Bootable USB stick](https://ubuntu.com/tutorials/install-ubuntu-desktop#3-create-a-bootable-usb-stick)
## Boot to BIOS setup screens
Press the F2 key on start up to enter the BIOS setup screens. Ensure that the BIOS is set to UEFI, and disable the Legacy option ROMS, and enable the secure boot. Here are a few refernce images on how to do it. Please be aware, your BIOS insterface can be a little different:

![image](https://user-images.githubusercontent.com/52679027/174858747-5383538b-bb34-4629-aed3-9b7e84392246.png)
![image](https://user-images.githubusercontent.com/52679027/174859187-b8cdaab2-f805-4841-b548-547c6f96fd42.png)
![image](https://user-images.githubusercontent.com/52679027/174859243-39573b11-b01e-43a2-a56b-aa2307b15b37.png)

After change BIOS, please save the change and exit

## Boot from USB flash drive
1. Insert the USB flash drice into the PC USB port and boot or restart the PC
2. Tap rapidly on teh F12 key when the logo appears during startup
3. Select USB device from the boot menu

## Install ubuntu 
Follow the [CANONICAL ubuntu step by step tutorials](https://ubuntu.com/tutorials/install-ubuntu-desktop#1-overview) to install Ubuntu 

## Formatting the second hard drive
A lot of time, you have one ssd drive, in which you intall your operating system, and another HDD drive to store your data. For the hard drive, you may want to delete the previous partitions, format it, mount it and then add it to the fstab. The following is a list of tutorials: 
* [How to Delete Partition in Linux](https://phoenixnap.com/kb/delete-partition-linux)
* [Format disk, mount the disk, update fstab](https://www.cyberciti.biz/faq/linux-disk-format/)

# Setup MinKNOW with Guppy GPU Basecaller
I was mainly following [MinKNOW tutrial](https://jhuapl-bio.github.io/Basestack/supplemental/minknow_guppy) to install MinKNOW, CUDA Toolkit, setup MinKNOW with Guppy GPU Basecaller

## Installing MiinKNOW

> wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
> echo "deb http://mirror.oxfordnanoportal.com/apt $(lsb_release -c | awk '{print $2}')-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
> sudo apt-get -y update
> sudo apt-get install -y minion-nc

## Install cuda
> wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
> sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
> wget https://developer.download.nvidia.com/compute/cuda/11.7.0/local_installers/cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
> sudo dpkg -i cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
> sudo cp /var/cuda-repo-ubuntu2004-11-7-local/cuda-*-keyring.gpg /usr/share/keyrings/
> sudo apt-get update
> sudo apt-get -y install cuda

reboot 
## Setup MinKNOW with Guppy GPU Baseller
> sudo mv /opt/ont/guppy/bin /opt/ont/guppy/bin.sav  &&    sudo mv /opt/ont/guppy/data /opt/ont/guppy/data.sav      # Save the old guppy just in case
> tar -xvzf ont-guppy_5.0.13_linux64.tar.gz 
> sudo cp -r ont-guppy/bin /opt/ont/guppy/bin && sudo cp -r ont-guppy/data /opt/ont/guppy/data # Move the newly downloaded guppy
> sudo /opt/ont/minknow/bin/config_editor --filename /opt/ont/minknow/conf/sys_conf --conf system --set on_acquisition_ping_failure=ignore
> sudo service minknow restart # Resart minknow


minion01_calgary@rice:~/Desktop$ sudo service guppyd stop
Warning: The unit file, source configuration file or drop-ins of guppyd.service changed on disk. Run 'systemctl daemon-reload' to reload units.


