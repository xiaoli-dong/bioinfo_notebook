
# How to Install Ubuntu Linux on your Dell Computer
https://www.dell.com/support/kbdoc/en-ca/000131655/how-to-install-ubuntu-linux-on-your-dell-pc
## Create a bootable USB flash drive
Download Ubuntu desktop image, the version we are using is [Ubuntu 22.04 LTS](https://ubuntu.com/download/desktop)
Follow the Canonical ubuntu totorials on how to [Create a Bootable USB stick](https://ubuntu.com/tutorials/install-ubuntu-desktop#3-create-a-bootable-usb-stick)
# Formatting the second hard drive
How to Delete Partition in Linux
https://phoenixnap.com/kb/delete-partition-linux

Format disk, mount the disk, update fstab
https://www.cyberciti.biz/faq/linux-disk-format/
sudo mkfs.ext4 /dev/sda

# Install Minknow
https://jhuapl-bio.github.io/Basestack/supplemental/minknow_guppy

# Install cuda
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.7.0/local_installers/cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo cp /var/cuda-repo-ubuntu2004-11-7-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda
reboot 

sudo mv /opt/ont/guppy/bin /opt/ont/guppy/bin.sav  &&    sudo mv /opt/ont/guppy/data /opt/ont/guppy/data.sav      # Save the old guppy just in case
tar -xvzf ont-guppy_5.0.13_linux64.tar.gz 
sudo cp -r ont-guppy/bin /opt/ont/guppy/bin && sudo cp -r ont-guppy/data /opt/ont/guppy/data # Move the newly downloaded guppy
#Disable online need for minknow to ping external servers
sudo /opt/ont/minknow/bin/config_editor --filename /opt/ont/minknow/conf/sys_conf --conf system --set on_acquisition_ping_failure=ignore
sudo service minknow restart # Resart minknow


minion01_calgary@rice:~/Desktop$ sudo service guppyd stop
Warning: The unit file, source configuration file or drop-ins of guppyd.service changed on disk. Run 'systemctl daemon-reload' to reload units.


