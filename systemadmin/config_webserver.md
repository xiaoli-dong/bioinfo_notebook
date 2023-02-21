
# Install and configure nginx web server
[reference](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/8/html/deploying_different_types_of_servers/setting-up-and-configuring-nginx_deploying-different-types-of-servers#installing-and-preparing-nginx_setting-up-and-configuring-nginx)
```
yum module list nginx  
yum install nginx  

firewall-cmd --permanent --add-port={80/tcp,443/tcp}  
firewall-cmd --reload  
systemctl enable nginx
systemctl start nginx
systemctl status nginx
```
#How do I turn SELinux off in Red Hat Enterprise Linux?  
https://access.redhat.com/solutions/3176  

[Install mongodb on redhat with SELinux enforced](https://technixleo.com/install-mongodb-on-centos-alma-rhel-9/)

chcon -v --type=httpd_sys_content_t /nfs/APL.../production/web  

# front-end
This application can only work with v16 

```
# if your system already have node npm installed, remove it first
yum remove -y nodejs npm
 curl -sL https://rpm.nodesource.com/setup_16.x | bash -

```

First, you need to add the NodeSource yum repository to your system. This can be done by using the following command as user root:
```
curl -sL https://rpm.nodesource.com/setup_16.x | bash -
 yum install nodejs
```
