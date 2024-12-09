
# Install Ubuntu on a Dell desktop configured for the Unified Extensible firmware interface (UEFI) BIOS

## Last updated on December 9, 2024

This is a tutorial on how to install Ubuntu on a dell PC, then configure MinKNOW, CUDA, and enable GPU based basecaller with MinKNOW. 

## Create a bootable USB flash drive
1. Download Ubuntu desktop image, the version we are using is [Ubuntu 20.04 LTS](https://ubuntu.com/download/desktop)
2. Follow the Canonical ubuntu totorials on how to [Create a Bootable USB stick](https://ubuntu.com/tutorials/install-ubuntu-desktop#3-create-a-bootable-usb-stick)
## Boot to BIOS setup screens
Press the F2 key on start up to enter the BIOS setup screens. Ensure that the BIOS is set to UEFI, and disable the Legacy option ROMS, and disable the secure boot. Here are a few reference images on how to do it. Please be aware, your BIOS insterface can be a little different:

![image](https://user-images.githubusercontent.com/52679027/174858747-5383538b-bb34-4629-aed3-9b7e84392246.png)
![image](https://user-images.githubusercontent.com/52679027/174859187-b8cdaab2-f805-4841-b548-547c6f96fd42.png)
![image](https://user-images.githubusercontent.com/52679027/174892849-f51e928b-cd52-48db-aa48-60b423083427.png)


After change BIOS, please save the change and exit

## Boot from USB flash drive
1. Insert the USB flash drice into the PC USB port and boot or restart the PC
2. Tap rapidly on the F12 key when the logo appears during startup
3. Select USB device from the boot menu

## Install ubuntu 
Follow the [CANONICAL ubuntu step by step tutorials](https://ubuntu.com/tutorials/install-ubuntu-desktop#1-overview) to install Ubuntu 

## Formatting the second hard drive
A lot of time, you have one ssd drive, in which you install your operating system, and another HDD drive to store your data. For the hard drive, you may want to delete the previous partitions, format it, mount it and then add it to the fstab. The following is a list of tutorials: 
* [How to Delete Partition in Linux](https://phoenixnap.com/kb/delete-partition-linux)
* [Format disk, mount the disk, update fstab](https://www.cyberciti.biz/faq/linux-disk-format/)
* [Disks utility to format a hard drive](https://www.wikihow.com/Format-a-Hard-Drive-Using-Ubuntu)

## Enabling SSH on Ubuntu
By default, when Ubuntu is first installed, remote access via SSH is not allowed. If you need it, enabling SSH on Ubuntu is fairly straightforward. Perform the following steps as root or user with sudo privileges to install and enable SSH on your Ubuntu system:

```
# Install openssh-server package
sudo apt update
sudo apt install openssh-server

# Verify that SSH is running by typing
sudo systemctl status ssh

# Ubuntu ships with a firewall configuration tool called UFW. If the firewall is enabled on your system, make sure to open the SSH port:
sudo ufw allow ssh

# connecting to the SSH Server
ssh username@ip_address
```

# Install CUDA Toolkit
Firstly, you need to ensure that your GPU is CUDA-capable by typing:

```
lspci | grep VGA
```
If you see your GPU model, for example: NVIDIA Corporation TU102 [GeForce RTX 2080 Ti] (rev A1) then you have a GPU available on your machine. IF you don’t see that and you know there is a GPU in the machine, try to install the drivers first. Then start to install CUDA Tookit, in my case, I installed [CUDA Toolkit 11.7](https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=Ubuntu&target_version=20.04&target_type=deb_local)

```
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.7.0/local_installers/cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo cp /var/cuda-repo-ubuntu2004-11-7-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda
```

Then, add these three lines to your $HOME/.bashrc

```
export LD_LIBRARY_PATH=/usr/local/cuda/lib64\
                         ${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export PATH=/usr/local/cuda/bin:$PATH
```

You should then **reboot your PC for cuda to take full effect**. Once rebooted, you should confirm that it is working by writing:
```
nvidia-smi
nvcc --version
```
* Notes: in my case, in the secure boot mode, those commands were not working properly. After dsiabling secure boot mode in BIOS, they worked.

Sample output from nvidia-smi command:  

![Screenshot from 2022-06-21 14-29-55](https://user-images.githubusercontent.com/52679027/174892221-0cda31c0-e1ea-4a88-abc7-27138a5ff8fe.png)

* Notes: before Guppy GPU Basecaller gets configured, the guppy_basecaller_server line will not shown in the screenshot

Sample output from nvcc --version command:  
![Screenshot from 2022-06-21 14-30-48](https://user-images.githubusercontent.com/52679027/174892486-3c303742-a0ff-4edd-b2ec-0056cdb9ed03.png)

## remove incorrectly installed nvidia and cuda toolkit
If you mess up the nvidia and CUDA ToolKit installation, here is some of the commands to remove nvidia, and cuda toolkits before a reinstallation.

```
 sudo apt-get --purge remove "*cublas*" "*cufft*" "*curand*"  "*cusolver*" "*cusparse*" "*npp*" "*nvjpeg*" "cuda*" "nsight*"
 sudo apt-get --purge remove "*nvidia*"
 sudo apt-get autoremove
 sudo reboot
```
After reboot, you can using command 
```
nvidia-smi
```
to make sure that your nvidia drive is working properly


# Install MinKNOW

For installation, you can just follow the offical document from nanopore. The new version of the MinKnow already includes the gpu based guppy. This has made the installation process much easier: 

```
wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://mirror.oxfordnanoportal.com/apt $(lsb_release -c | awk '{print $2}')-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get -y update
sudo apt-get install -y minion-nc
```

## Edit app_conf file  
Before editing, I always backup the original app_conf file to app_conf.bak. You will also need to modifying /opt/ont/minknow/conf/app_conf file. Adjust the "gpu_calling" field from false to true, being careful not to modify/delete any commas or quotations. See the reference image below:  

![Screenshot from 2022-06-21 14-02-38](https://user-images.githubusercontent.com/52679027/174887939-b7b24cfd-54f0-4191-bc7c-791fa767f3ce.png)

For the projects we are doing, we are very often required to check barcodes at both ends during the barcode resolving stage. I used config_editor added that requirement to the minknow configurtion file app_conf

```
sudo /opt/ont/minknow/bin/config_editor --conf application --filename /opt/ont/minknow/conf/app_conf --set guppy.server_config.extra_arguments="--require_barcodes_both_ends"

```

## Customize the data output directory
By default, MinKNOW will output the data to /var/lib/minknow/data directory and output the log files to /var/log/minknow directory. If you want to change the default data and log output directory, you will need to edit /opt/ont/minknow/conf/user_conf file. Before editing, please backup the original file using the command below:

```
sudo cp -v user_conf user_conf.bak

```
You can refer to the screenshots on how to modify the user_conf file

<figure>
  <figcaption>The ouput_dirs section of the original user_conf</figcaption>
  <img
  src="https://user-images.githubusercontent.com/52679027/186774882-638fc3ff-9703-422e-b9c9-d56bd2f40bcc.png"
  alt="The ouput_dirs section of the original user_conf">
  
</figure>  



<figure>
  <figcaption>The ouput_dirs section of the user_conf after customization</figcaption>
  <img
  src="https://user-images.githubusercontent.com/52679027/186774884-c5e83cba-71ad-439c-a638-641bb1351a6c.png"
  alt="The ouput_dirs section of the user_conf after customization">
  
</figure>

After modifying the user_conf file, restart MinKNOW, Guppyd

```
sudo service minknow stop # Resart minknow
sudo service guppyd stop
sudo service guppyd start
sudo service minknow restart
```

When I connected the device to the computer, I got the error message in the log file as " ERROR: libusb: error [get_usbfs_fd] libusb couldn't open USB device". See the screenshot:

![Screenshot 2022-11-25 160147](https://user-images.githubusercontent.com/52679027/204061870-62ba26c6-01df-4c49-b2f6-a1824c0fa87f.png)

Then, I changed "User=minkow Group=minknow" to "User=root Group=root". After restart the service, then the error message has gone!!!!!!!!!!!!!

# Troubleshooting reference
(GPU Calling in MinKNOW)[https://gringer.gitlab.io/presentation-notes/2021/10/08/gpu-calling-in-minknow/]

## Upgrade to the newer verion of minknow
### upgrade Ubuntu
Fully update the system. The upgrade process works best when the current system has all the latest updates installed. You should confirm that these commands complete successfully and that no further updates are available. We also suggest rebooting the system after all the updates are applied, to ensure the latest kernel is being run. To upgrade, run the following commands:

```
apt-get auto-remove && apt-get clean && apt-get update && apt-get upgrade do-release-upgrade
```
### install CUDA Toolkit and driver 12.6

#### clean the previous installation
```
apt-get --purge remove -y "*cublas*" "*cufft*" "*curand*"  "*cusolver*" "*cusparse*" "*npp*" "*nvjpeg*" "cuda*" "nsight*"
apt-get --purge remove -y "*nvidia*"
apt-get autoremove -y
reboot
```

#### install CUDA Toolkit 12.6
```
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/12.6.3/local_installers/cuda-repo-ubuntu2004-12-6-local_12.6.3-560.35.05-1_amd64.deb

cp /var/cuda-repo-ubuntu2004-12-6-local/cuda-F2089CC5-keyring.gpg /usr/share/keyrings/
dpkg -i cuda-repo-ubuntu2004-12-6-local_12.6.3-560.35.05-1_amd64.deb
cp /var/cuda-repo-ubuntu2004-12-6-local/cuda-*-keyring.gpg /usr/share/keyrings/
apt-get update
#apt-get -y install cuda-toolkit-12-6
apt-get -y install cuda
reboot
```
You should then **reboot your PC for cuda to take full effect**. Once rebooted, you should confirm that it is working by writing:

```
nvidia-smi
nvcc --version
```

### Install MinKnow Version 24.06.16
This version of the minknow integrated Dorado into the MinKNOW

#### remove the previouse installation
```
apt-get purge -y ont-*
apt-get autoremove
```

#### For Ubuntu 20 to add the Oxford Nanopore apt repository, run the command below on a terminal window:
```
sudo apt update
sudo apt install wget
wget -O- https://cdn.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://cdn.oxfordnanoportal.com/apt focal-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
```
Install GPU version of the MinKNOW using the command:
```
sudo apt update
sudo apt install ont-standalone-minknow-gpu-release
reboot
```
## References
[GPU Calling in MinKNOW](https://gringer.gitlab.io/presentation-notes/2021/10/08/gpu-calling-in-minknow/)

[How to Install Ubuntu Linux on your Dell Computer](https://www.dell.com/support/kbdoc/en-ca/000131655/how-to-install-ubuntu-linux-on-your-dell-pc)

[Ubuntu Wiki releases page](https://wiki.ubuntu.com/Releases?_ga=2.181102951.1543743502.1714578127-482619844.1714578127)
