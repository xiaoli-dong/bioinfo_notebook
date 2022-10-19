# Hardware optimization for the slurmctld master server
[SchedMD](https://www.schedmd.com/) recommends that the slurmctld server should have only a few, but very fast CPU cores, in order to ensure the best responsiveness. The file system for /var/spool/slurmctld/ should be mounted on the fastest possible disks (SSD or NVMe if possible).
# Create global user accounts
There must be a uniform user and group name space (including UIDs and GIDs) across the cluster. Slurm and MUNGE require consistent UID and GID across all servers and nodes in the cluster, including the slurm and munge user and make sure that these same users are created identically on all nodes. This must be done prior to installing RPMs (which would create random UID/GID pairs if these users don’t exist).

# MUNGE authentication service
* Install the MUNGE RPM packages 
# MUNGE configuration and testing
* On the Head/Master node (only) create a secret key to be used globally on every node
* Securely propagate /etc/munge/munge.key (e.g., via SSH) to all other hosts within the same security realm:
* Make sure to set the correct ownership and mode on all nodes:
# Then enable and start the MUNGE service on all nodes:



https://wiki.fysik.dtu.dk/Niflheim_system/Slurm_installation/
