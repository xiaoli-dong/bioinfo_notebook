# How to install and configure R+Rstudio server on RHEL 8 Linux System

## Install R from source code

Following the [Posit Documentation](https://docs.posit.co/resources/install-r-source/): 
```
# Enable the Extra Packages for Enterprise Linux (EPEL) repository
yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm

# Enable the CodeReady Linux Builder repository:
subscription-manager repos --enable codeready-builder-for-rhel-8-x86_64-rpms

# install the build dependencies for R:
dnf builddep R

# download and extract R
export R_VERSION=4.3.2
curl -O https://cran.rstudio.com/src/base/R-4/R-${R_VERSION}.tar.gz
tar -xzvf R-${R_VERSION}.tar.gz
cd R-${R_VERSION}

# build and install R
#--enable-R-shlib	Required to use R with RStudio.
#--enable-memory-profiling	Enables support for Rprofmem() and tracemem(), used to measure memory use in R code.
 ./configure --prefix=/nfs/APL_Genomics/apps/production/R/R-4.3.2/build --enable-R-shlib --enable-memory-profiling
make -j 8
make install
```

If you need to install R packages on an individual basis please create a directory within your preferred folder and then use a variable to point R to the directory, so the package can be installed. At the shell prompt, enter the following command:

```
export RHOME="${prog}/R"
#Pre-create a directory (such as ~/R-packages) and set the R_LIBS environment variable before installing R packages
export R_LIBS="${prog}/R/R-packages"
```

## Install and configure RStudio server
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

However, we cannot leave our SELinux in the permissive mode. Referencing the [forum post](https://github.com/rstudio/rstudio/issues/4937) here: it seems to indicate that the issue is with the selinux context on the rstudio binaries themselves. Typically binaries are installed in the /usr/bin, /usr/sbin, etc., and not the/usr/lib directory, or sub-directories.  By default files/directories created under the /usr/lib directory have a context of lib_t: e.g. from /usr/lib/

```
# ls -lZ | grep rstudio
drwxr-xr-x.  9 root root system_u:object_r:lib_t:s0              190 Feb  1 08:44 rstudio-server
```
whereas the expected security context of binaries would be bin_t: e.g. from /usr/bin
```
# ls -lZ | grep zless
lrwxrwxrwx. 1 root root    system_u:object_r:bin_t:s0                            6 Aug 12  2018 bzless -> bzmore
-rwxr-xr-x. 1 root root    system_u:object_r:bin_t:s0                            1802 May 31  2022 xzless
-rwxr-xr-x. 1 root root    system_u:object_r:bin_t:s0                            2205 Apr 20  2022 zless
```
The forum post above seems to be pretty vocal about how this is an issue with the rstudio product (and more how it’s packaged/installed) and is something the developers of rstudio need to address to properly accommodate selinux enabled distributions, it’s not something you could have addressed or updated via a config file as part of your installation.
As suggested in the article we did the following:
```
# Changed the selinux context of the binary files using the semanage command:
$ semanage fcontext -a -t bin_t '/usr/lib/rstudio-server/bin(/.*)?'
Restored the permissions using the restorecon command
$ restorecon -r /usr/lib/rstudio-server/bin/
Restarted the rstudio server service
$ systemctl restart rstudio-server.service
$ systemctl status rstudio-server.service
Put SELinux back into enforcing mode
$ getenforce
$ setenforce 1
$ sestatus
```
Had users login to confirm things were “working” and checked the journald logs for selinux issues
```
$ journalctl -t setroubleshoot
```
We searched for any rstudio or rsession lines and none were found after we updated the context to bin_t for the rstudio binaries
After these steps we were able to confirm that our users were able to login to the rstudio UI/webpage, after first logging into the server via ssh so their home directory would be created.



