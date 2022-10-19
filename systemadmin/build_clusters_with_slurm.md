# Hardware optimization for the slurmctld master server
[SchedMD](https://www.schedmd.com/) recommends that the slurmctld server should have only a few, but very fast CPU cores, in order to ensure the best responsiveness. The file system for /var/spool/slurmctld/ should be mounted on the fastest possible disks (SSD or NVMe if possible).
