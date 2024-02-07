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
# cat <<-EOF >> /etc/profile.d/R.sh
export PATH=/nfs/APL_Genomics/apps/production/R/bin:\$PATH
EOF


```

https://docs.posit.co/resources/install-r/

https://docs.posit.co/resources/install-r-source/

## install RStudio Server
https://web.mit.edu/r/current/RStudio/INSTALL

https://higgi13425.github.io/medical_r/posts/2020-12-06-setting-up-a-rstudio-server-with-free-software-version/

Download RStudio sever: https://posit.co/download/rstudio-server/

```
wget wget https://download2.rstudio.org/server/rhel8/x86_64/rstudio-server-rhel-2023.12.1-402-x86_64.rpm
yum install rstudio-server-rhel-2023.12.1-402-x86_64.rpm

```
You need to edit the configuration file /etc/rstudio/rserver.conf to include the following line.

As we did earlier, make sure this path is the same path where you installed R.

```
rsession-which-r=/nfs/APL_Genomics/apps/production/R/bin/R

if R directory is symbolink
ERROR Unable to determine real path of R script /nfs/APL_Genomics/apps/production/R/bin/R (system error 13 (Permission denied));
if point to the realpaht
ERROR Error reading R script (/nfs/APL_Genomics/apps/production/R/R-4.3.2/build/bin/R), system error 2 (No such file or directory)

After copyt the content in the build directory to /opt/R, it started properly

sudo systemctl daemon-reload 
sudo systemctl start rstudio-server 
sudo systemctl enable rstudio-server

 firewall-cmd --permanent --zone=public --add-port=8787/tcp
firewall-cmd --reload
http://10.106.109.188:8787/
```

![image](https://github.com/xiaoli-dong/bioinfo_notebook/assets/52679027/3f80250a-e93a-4a90-9c82-f7d0d283e0c6)

When login to the RStudio server with the LDAP username and password, we get the Following error message:

"Error: Incorrect or invalid username/password"

RStudio connects to LDAP via PAM and I used pamtester to identify login issues: 
```
pamtester --verbose rstudio <username> authenticate
pamtester: authentication failed
```



To integrate RStudio server with PAM. I copied the PAM profile to use with RStudio. Reference link is here [Using LDAP authentication with RStudio Workbench / RStudio Server Pro](https://support.posit.co/hc/en-us/articles/232226708-Using-LDAP-authentication-with-RStudio-Workbench-RStudio-Server-Pro)

```
cp /etc/pam.d/login /etc/pam.d/rstudio

pamtester --verbose rstudio my_userid authenticate acct_mgmt setcred open_session close_session
pamtester: invoking pam_start(rstudio, my_userid, ...)
pamtester: performing operation - authenticate
Password:
pamtester: successfully authenticated
pamtester: performing operation - acct_mgmt
pamtester: account management done.
pamtester: performing operation - setcred
pamtester: credential info has successfully been set.
pamtester: performing operation - open_session
pamtester: sucessfully opened a session
pamtester: performing operation - close_session
pamtester: session has successfully been closed.
```
Tried to access RStudio server again with: http://10.106.109.188:8787/ using LDAP username and password, get the same error message:
"Error: Incorrect or invalid username/password"

Because our system is using SELinux, I temporarily set it to permissive according to the procedures here: [Using LDAP authentication with RStudio Workbench / RStudio Server Pro](https://support.posit.co/hc/en-us/articles/15173704481943-Active-Directory-LDAP-user-not-able-to-login-permission-denied-on-PAM-acct-mgmt)

```
setenforce 0
systemctl restart rstudio-server 
```
Now, I can access the RStudio server

