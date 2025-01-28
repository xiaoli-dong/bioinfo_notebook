
# How to setup a HPC cluster using slurm
Last updated:  January 28, 2025

## Hardware optimization for the slurmctld master server
[SchedMD](https://www.schedmd.com/) recommends that the slurmctld server should have only a few, but very fast CPU cores, in order to ensure the best responsiveness. The file system for /var/spool/slurmctld/ should be mounted on the fastest possible disks (SSD or NVMe if possible).
## Create global user accounts
There must be a uniform user and group name space (including UIDs and GIDs) across the cluster. Slurm and MUNGE require consistent UID and GID across all servers and nodes in the cluster, including the slurm and munge user and **make sure that these same users are created identically on all nodes**. This must be done prior to installing RPMs (which would create random UID/GID pairs if these users don’t exist).

It is very important to avoid UID and GID below 1000, as defined in the standard configuration file /etc/login.defs by the parameters UID_MIN, UID_MAX, GID_MIN, GID_MAX, see also https://en.wikipedia.org/wiki/User_identifier.

```
# run the commands on all the nodes (head and computer nodes)
export MUNGEUSER=1005
groupadd -g $MUNGEUSER munge
useradd  -m -c "MUNGE Uid 'N' Gid Emporium" -d /var/lib/munge -u $MUNGEUSER -g munge  -s /sbin/nologin munge
export SlurmUSER=1001
groupadd -g $SlurmUSER slurm
useradd  -m -c "Slurm workload manager" -d /var/lib/slurm -u $SlurmUSER -g slurm  -s /bin/bash slurm
```

## MUNGE packages installation and configuration

### What is Munge?
Munge is a service that uses a shared secret to generate cryptographic tokens (authentication tokens) for identifying users and ensuring trust across different nodes. It encrypts and signs the data to prevent unauthorized access and impersonation of users.

### Role of Munge in Slurm
In a Slurm-based cluster: 
* Munge provides authentication for communication between different Slurm daemons across nodes in the cluster.
* Munge ensures that job submissions and commands (e.g., sbatch, scontrol, srun) are securely authenticated.
* It prevents impersonation attacks, where an unauthorized user might try to execute jobs or control Slurm daemons.
* Slurm relies on Munge to verify the identity of the user submitting jobs or managing the cluster.
### Install MUNGE
To use Munge with Slurm, the Munge service must be installed, configured, and running on all nodes in the cluster. Here's a step-by-step guide to set up Munge in a Slurm cluster:

```
#PowerTools is a CentOS repository. On RHEL 8 we have the CodeReady
subscription-manager repos --enable codeready-builder-for-rhel-8-x86_64-rpms
yum install munge munge-libs munge-devel

```

### MUNGE configuration and testing

* On the Head/Master node (only) create a secret key to be used globally on every node
* Securely propagate /etc/munge/munge.key (e.g., via SSH) to all other hosts within the same security realm:
* Make sure to set the correct ownership and mode on all nodes:

#### Create MUNGE secret key
```
sudo yum install rng-tools -y
sudo rngd -r /dev/urandom 
sudo /usr/sbin/create-munge-key -r -f

# This would be 1024 bytes of random data printed to munge.key file
sudo sh -c  "dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key"
sudo chown munge: /etc/munge/munge.key

# 400 only ower can read the content
# This ensures that the private key is readable only by the owner, protecting it from being accessed by others.
sudo chmod 400 /etc/munge/munge.key
```
#### Copy all munge.key to all the nodes and make sure to set the correct ownership and mode on all nodes:

```
cp -p /your_path_to/munge.key /etc/munge/
mkdir /var/log/munge
chown munge: /etc/munge/munge.key
chmod 400 /etc/munge/munge.key
chown -R munge: /etc/munge/ /var/log/munge/

#Grant full access (read, write, and execute) to the owner of the file or directory
#Deny all access to the group and others.
chmod 0700 /etc/munge/ /var/log/munge/
```
### Enable MUNGE authentication service on all the nodes

```
sudo systemctl enable munge
sudo systemctl start munge
```


## Build slurm on the head node

### Enable slurm accounting database 

```
sudo yum install mariadb-server mariadb-devel -y
```

### Install Slurm prerequisites as well as several optional packages that enable Slurm plugins as described in the Slurm_Quick_Start guide:
```
yum install rpm-build gcc python3 openssl openssl-devel pam-devel numactl numactl-devel hwloc hwloc-devel munge munge-libs munge-devel lua lua-devel readline-devel rrdtool-devel ncurses-devel gtk2-devel libibmad libibumad perl-Switch perl-ExtUtils-MakeMaker xorg-x11-xauth http-parser-devel json-c-devel

# this is required by libssh2-devel man2html
yum module enable virt-devel

yum install libssh2-devel man2html

#If you want to build the Slurm REST API daemon named slurmrestd (from Slurm 20.02 and newer), or if you want to use the slurm.conf ResumeProgram and SuspendProgram from the Power_Saving_Guide, then you make sure to install these prerequisites before building RPMs:
yum install http-parser-devel json-c-devel
```

### Build yum rpm package

```
mkdir slurm-tmp
cd slurm-tmp
export VER=22.05.5
wget https://download.schedmd.com/slurm/slurm-$VER.tar.bz2
rpmbuild -ta slurm-$VER.tar.bz2  
rm slurm-$VER.tar.bz2
cd ..
rmdir slurm-tmp 
```

## install and configure slurm on head node

```
# get perl-Switch
yum install cpan -y 
cd ~/rpmbuild/RPMS/x86_64/
yum --nogpgcheck localinstall slurm-22.05.5-1.el8.x86_64.rpm slurm-contribs-22.05.5-1.el8.x86_64.rpm slurm-devel-22.05.5-1.el8.x86_64.rpm slurm-example-configs-22.05.5-1.el8.x86_64.rpm  slurm-libpmi-22.05.5-1.el8.x86_64.rpm  slurm-openlava-22.05.5-1.el8.x86_64.rpm slurm-pam_slurm-22.05.5-1.el8.x86_64.rpm  slurm-perlapi-22.05.5-1.el8.x86_64.rpm  slurm-slurmctld-22.05.5-1.el8.x86_64.rpm  slurm-slurmd-22.05.5-1.el8.x86_64.rpm slurm-slurmdbd-22.05.5-1.el8.x86_64.rpm slurm-torque-22.05.5-1.el8.x86_64.rpm -y

## configure slurm
sudo mkdir /var/spool/slurm
sudo chown slurm:slurm /var/spool/slurm
sudo chmod 755 /var/spool/slurm
sudo mkdir /var/spool/slurm/slurmctld
sudo chown slurm:slurm /var/spool/slurm/slurmctld
sudo chmod 755 /var/spool/slurm/slurmctld
sudo mkdir -p /var/spool/slurm/cluster_state
sudo chown slurm:slurm /var/spool/slurm/cluster_state```
 
 #on login node
 firewall-cmd --permanent --zone=public --add-port=6817/udp
 firewall-cmd --permanent --zone=public --add-port=6817/tcp
 firewall-cmd --permanent --zone=public --add-port=6818/tcp
 firewall-cmd --permanent --zone=public --add-port=6818/udp
 firewall-cmd --permanent --zone=public --add-port=7321/udp
 firewall-cmd --permanent --zone=public --add-port=7321/tcp
 firewall-cmd --permanent --zone=public --add-port=6819/tcp
 firewall-cmd --permanent --zone=public --add-port=6819/udp
 firewall-cmd --add-port=30000-60000/udp --permanent
 firewall-cmd --add-port=30000-60000/tcp --permanent

sudo touch /var/log/slurmctld.log
sudo chown slurm:slurm /var/log/slurmctld.log
sudo touch /var/log/slurm_jobacct.log /var/log/slurm_jobcomp.log
sudo chown slurm: /var/log/slurm_jobacct.log /var/log/slurm_jobcomp.log


#we are doing Configless Slurm setup by adding "SlurmctldParameters=enable_configless" to the slurm.conf
#https://slurm.schedmd.com/configless_slurm.html


# enable and start services
systemctl enable slurmctld
systemctl enable slurmdbd
systemctl start slurmctld.service
```
### Update the slurm.conf file
When the slurm.conf file is updated, please run the following command to instruct all Slurm daemons to re-read the configuration file
```
scontrol reconfigure
```
Otherwise, you will get the below error message
```
error: Node xxx appears to have a different slurm.conf than the slurmctld.
```
## Install and configure slurm on all the computing nodes
```
# make a copy of the slurm rpm build from the head node before start
cd ~/rpmbuild/RPMS/x86_64/

# isntall slurm
yum --nogpgcheck localinstall slurm-22.05.5-1.el8.x86_64.rpm slurm-contribs-22.05.5-1.el8.x86_64.rpm slurm-devel-22.05.5-1.el8.x86_64.rpm slurm-example-configs-22.05.5-1.el8.x86_64.rpm  slurm-libpmi-22.05.5-1.el8.x86_64.rpm  slurm-openlava-22.05.5-1.el8.x86_64.rpm slurm-pam_slurm-22.05.5-1.el8.x86_64.rpm  slurm-perlapi-22.05.5-1.el8.x86_64.rpm  slurm-slurmctld-22.05.5-1.el8.x86_64.rpm  slurm-slurmd-22.05.5-1.el8.x86_64.rpm slurm-slurmdbd-22.05.5-1.el8.x86_64.rpm slurm-torque-22.05.5-1.el8.x86_64.rpm -y

# configure slurm
mkdir -p /var/spool/slurm/slurmd
mkdir -p /var/log/slurm
touch /var/log/slurm/slurmd.log
chown -R slurm:slurm /var/spool/slurm /var/log/slurm
mkdir -p /run/slurm
chown -R slurm:slurm /run/slurm
```
firewall will block connections between nodes so in case of cluster with multiple nodes adapt the firewall on the compute nodes, we will disable the firewall on the computing node

```
systemctl stop firewalld.service
systemctl disable firewalld.service
```
otherwise, your computing nodes will be in "down" state after it stay in idle for a while

![image](https://user-images.githubusercontent.com/52679027/200654143-e745ee97-eacd-4f46-849c-1dd3fabfa5dc.png)

When using "enable_configles" option, you must configure the slurmd to get its configs from the slurmctld. This can be accomplished by launching slurmd with the "--conf-server" option in the slurmd.service file:
 ```
 ExecStart=/usr/sbin/slurmd --conf-server your_head_node_hostname:6817 -D -s
 ```
 Then you can enbale and start the slurmd service
 ```
 systemctl daemon-reload,
 systemctl enable slurmd.service
 systemctl start slurmd.service
 systemctl status slurmd.service
 ```
    
 make sure to restart munge.service


## Account setup
```
sacctmgr add User Accounts=all xiaolidong
sacctmgr add account all Description="all" Organization=all
sacctmgr show User
sacctmgr show association
sacctmgr list account
```
## Restart a node
After a node shutdown and upgrade, the node went to state down. The following commands will bring the node from down state to idle state

```
#If no jobs are currently running on the node:
scontrol update nodename=computer_node_name state=idle

#If jobs are running on the node:
scontrol update nodename=computer_node_name state=resume
```


## References

* [Tasks for Account Coordinators](https://rcic.uci.edu/hpc3/account-control.html)
* [Slurm administration](https://documentation.tjhsst.edu/services/cluster/slurm-administration)
* [faq](https://hpc.pku.edu.cn/_book/guide/faq.html)
* [Slurm installation and upgrading](https://wiki.fysik.dtu.dk/Niflheim_system/Slurm_installation/)
* [AUTOMATIC SLURM BUILD AND INSTALLATION SCRIPT](https://www.ni-sp.com/slurm-build-script-and-container-commercial-support/)
