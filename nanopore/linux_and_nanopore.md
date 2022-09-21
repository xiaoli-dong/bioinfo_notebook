This is a tutorial on how to install Ubuntu on a dell PC, then configure MinKNOW, CUDA, and enable GPU based guppy basecaller with MinKNOW  
# Install Ubuntu on a Dell desktop configured for the Unified Extensible firmware interface (UEFI) BIOS
I was using [How to Install Ubuntu Linux on your Dell Computer](https://www.dell.com/support/kbdoc/en-ca/000131655/how-to-install-ubuntu-linux-on-your-dell-pc) as a reference for Ubuntu setup

## Create a bootable USB flash drive
1. Download Ubuntu desktop image, the version we are using is [Ubuntu 22.04 LTS](https://ubuntu.com/download/desktop)
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

# Setup MinKNOW with Guppy GPU Basecaller
I was mainly using [MinKNOW tutorial](https://jhuapl-bio.github.io/Basestack/supplemental/minknow_guppy) as a reference to install MinKNOW, CUDA Toolkit, setup MinKNOW with Guppy GPU Basecaller

## Install MinKNOW
```
wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://mirror.oxfordnanoportal.com/apt $(lsb_release -c | awk '{print $2}')-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get -y update
sudo apt-get install -y minion-nc
```
## Install cuda
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

Sample output from nvcc --version command:  
![Screenshot from 2022-06-21 14-30-48](https://user-images.githubusercontent.com/52679027/174892486-3c303742-a0ff-4edd-b2ec-0056cdb9ed03.png)


## Setup MinKNOW with Guppy GPU Baseller
Finally, we can start to configure MinKNOW to use a GPU-capable version of guppy and to make the guppy basecaller plays nice with the installed MinKNOW.

```
/opt/ont/minknow/guppy/bin/guppy_basecaller --version
```
You should see a version, for example, we are using 6.1.5. You MUST download the same version of the Guppy:

```
# Save the old guppy just in case
sudo mv /opt/ont/guppy/bin /opt/ont/guppy/bin.sav && sudo mv /opt/ont/guppy/data /opt/ont/guppy/data.sav

tar -xvzf ont-guppy_6.1.5_linux64.tar.gz

# Move the newly downloaded guppy
sudo cp -r ont-guppy/bin /opt/ont/guppy/bin && sudo cp -r ont-guppy/data /opt/ont/guppy/data

#Disable online need for MinKNOW to ping external servers
sudo /opt/ont/minknow/bin/config_editor --filename /opt/ont/minknow/conf/sys_conf --conf system --set on_acquisition_ping_failure=ignore

# restart MinKNOW
sudo service minknow restart
```

## Edit the existing guppyd service file using 
```
sudo vi /etc/systemd/system/guppyd.service

# Change the line
> ExecStart=/opt/ont/guppy/bin/guppy_basecall_server --log_path /var/log/guppy --config dna_r9.4.1_450bps_fast.cfg --num_callers 1 --cpu_threads_per_caller 2 --port /tmp/.guppy/5555 --ipc_threads 3
# to 
> ExecStart=/opt/ont/guppy/bin/guppy_basecall_server --log_path /var/log/guppy --config dna_r9.4.1_450bps_fast.cfg --port /tmp/.guppy/5555 -x cuda:all
```

## Edit the existing override.conf file
```
sudo vi /etc/systemd/system/guppyd.service.d/override.conf

# Change the line
> ExecStart=/opt/ont/guppy/bin/guppy_basecall_server --log_path /var/log/guppy --config dna_r9.4.1_450bps_fast.cfg --num_callers 1 --cpu_threads_per_caller 2 --port /tmp/.guppy/5555 --ipc_threads 3
# to 
> ExecStart=/opt/ont/guppy/bin/guppy_basecall_server --log_path /var/log/guppy --config dna_r9.4.1_450bps_fast.cfg --port /tmp/.guppy/5555 -x cuda:all
```

## Edit app_conf file  
Before editing, I always backup the original app_conf file to app_conf.backup. You will also need to modifying /opt/ont/minknow/conf/app_conf file. Adjust the "gpu_calling" field from false to true, being careful not to modify/delete any commas or quotations. See the reference image below:  

![Screenshot from 2022-06-21 14-02-38](https://user-images.githubusercontent.com/52679027/174887939-b7b24cfd-54f0-4191-bc7c-791fa767f3ce.png)

For the projects we are doing, we are very often required to check barcodes at both ends during the barcode resolving stage. I used config_editor added that requirement to the minknow configurtion file app_conf

```
sudo /opt/ont/minknow/bin/config_editor --conf application --filename /opt/ont/minknow/conf/app_conf --set guppy.server_config.extra_arguments="--require_barcodes_both_ends"

```

After all the changes, restart MinKNOW, Guppyd

```
sudo service minknow stop # Resart minknow
sudo service guppyd stop
sudo service guppyd start
sudo service minknow restart
```
From there you are all set to run basecalling directly within the MinKNOW application.

## Customize the data output directory
By default, MinKNOW will output the data to /var/lib/minknow/data directory and output the log files to /var/log/minknow directory. If you want to change the default data and log output directory, you will need to edit /opt/ont/minknow/conf/user_conf file. Before editing, please backup the original file using the command below:

```
sudo cp -v user_conf user_conf.backup 

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
