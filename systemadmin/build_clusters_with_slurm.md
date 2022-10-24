# Hardware optimization for the slurmctld master server
[SchedMD](https://www.schedmd.com/) recommends that the slurmctld server should have only a few, but very fast CPU cores, in order to ensure the best responsiveness. The file system for /var/spool/slurmctld/ should be mounted on the fastest possible disks (SSD or NVMe if possible).
# Create global user accounts
There must be a uniform user and group name space (including UIDs and GIDs) across the cluster. Slurm and MUNGE require consistent UID and GID across all servers and nodes in the cluster, including the slurm and munge user and make sure that these same users are created identically on all nodes. This must be done prior to installing RPMs (which would create random UID/GID pairs if these users don’t exist).

It is very important to avoid UID and GID below 1000, as defined in the standard configuration file /etc/login.defs by the parameters UID_MIN, UID_MAX, GID_MIN, GID_MAX, see also https://en.wikipedia.org/wiki/User_identifier.

```
export MUNGEUSER=1005
groupadd -g $MUNGEUSER munge
useradd  -m -c "MUNGE Uid 'N' Gid Emporium" -d /var/lib/munge -u $MUNGEUSER -g munge  -s /sbin/nologin munge
export SlurmUSER=1001
groupadd -g $SlurmUSER slurm
useradd  -m -c "Slurm workload manager" -d /var/lib/slurm -u $SlurmUSER -g slurm  -s /bin/bash slurm
```

# install munge packages
```
#PowerTools is a CentOS repository. On RHEL 8 we have the CodeReady
subscription-manager repos --enable codeready-builder-for-rhel-8-x86_64-rpms
yum install munge munge-libs munge-devel

```
# MUNGE configuration and testing

* On the Head/Master node (only) create a secret key to be used globally on every node
* Securely propagate /etc/munge/munge.key (e.g., via SSH) to all other hosts within the same security realm:
* Make sure to set the correct ownership and mode on all nodes:

```
sudo yum install rng-tools -y
sudo rngd -r /dev/urandom
 
sudo /usr/sbin/create-munge-key -r -f
 
sudo sh -c  "dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key"
sudo chown munge: /etc/munge/munge.key
sudo chmod 400 /etc/munge/munge.key
```
copy all munge.key to all the nodes and make sure to set the correct ownership and mode on all nodes:
```
cp -p /nfs/APL_Genomics/munge.key /etc/munge/
mkdir /var/log/munge
chown munge: /etc/munge/munge.key
chmod 400 /etc/munge/munge.key
chown -R munge: /etc/munge/ /var/log/munge/
chmod 0700 /etc/munge/ /var/log/munge/
```
Enable MUNGE authentication service on all the nodes
```
sudo systemctl enable munge
sudo systemctl start munge
```



https://wiki.fysik.dtu.dk/Niflheim_system/Slurm_installation/
