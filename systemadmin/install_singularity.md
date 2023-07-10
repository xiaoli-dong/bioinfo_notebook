```
yum update -y && \
    yum groupinstall -y 'Development Tools' && \
    yum install -y \
    openssl-devel \
    libuuid-devel \
    libseccomp-devel \
    wget \
    squashfs-tools

 wget https://go.dev/dl/go1.20.5.linux-amd64.tar.gz
tar -C /nfs/APL_Genomics/apps/production/go -xvzf go1.20.5.linux-amd64.tar.gz
```
Download and install Singularity from a release
```
wget https://github.com/sylabs/singularity/releases/download/v3.11.4/singularity-ce-3.11.4.tar.gz
tar -xvzf singularity-ce-3.11.4.tar.gz
cd singularity-ce-3.11.4/
./mconfig  --without-suid --prefix=/nfs/APL_Genomics/apps/production/singularity/
make -C ./builddir
make -C ./builddir install
```
modify your bashrc file
```
### configularity singularity to be shared across all the nodes
export SINGULARITY="${prog}/singularity"
export SINGULARITY_HOME="${prog}/singularity_home"
export SINGULARITY_CACHEDIR="${SINGULARITY_HOME}/SINGULARITY_CACHEDIR"
export SINGULARITY_TMPDIR="${SINGULARITY_HOME}/SINGULARITY_TMPDIR"
export SINGULARITY_PULLFOLDER="${SINGULARITY_HOME}/SINGULARITY_PULLFOLDER"
export SINGULARITY_LOCALCACHEDIR="${SINGULARITY_HOME}/SINGULARITY_LOCALCACHEDIR"
export SINGULARITY_BIND="/nfs/APL_Genomics,/nfs/Genomics_DEV"

```
