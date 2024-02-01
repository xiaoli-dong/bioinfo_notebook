# How to install and configure R and Rstudio server on RHEL 8 Linux System

## Install R from source code 
```
# Enable the Extra Packages for Enterprise Linux (EPEL) repository
yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm
#Enable the CodeReady Linux Builder repository:
sudo subscription-manager repos --enable codeready-builder-for-rhel-8-x86_64-rpms
# install the build dependencies for R:
dnf builddep R

export R_VERSION=4.3.2
curl -O https://cran.rstudio.com/src/base/R-4/R-${R_VERSION}.tar.gz
tar -xzvf R-${R_VERSION}.tar.gz
cd R-${R_VERSION}

#--enable-R-shlib	Required to use R with RStudio.
#--enable-memory-profiling	Enables support for Rprofmem() and tracemem(), used to measure memory use in R code.
 ./configure --prefix=/nfs/APL_Genomics/apps/production/R/R-4.3.2/build --enable-R-shlib --enable-memory-profiling
make -j 8

export RHOME="${prog}/R"
#Pre-create a directory (such as ~/R-packages) and set the R_LIBS environment variable before installing R packages. (See sample steps below.)
export R_LIBS="${prog}/R/R-packages"

# add R to the path

```

https://docs.posit.co/resources/install-r/
https://docs.posit.co/resources/install-r-source/
