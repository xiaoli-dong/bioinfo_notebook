On systems running SELinux, all processes and files are labeled in a way that represents security-relevant information. This information is called the SELinux context. For files, this is viewed using the 
```
ls -Z
```

 If you want to allow nginx to have append access on the nginx-error.log file
                                                        Then you need to change the label on '/apps/web/logs/nginx/nginx-error.log'
                                                        Do
                                                        # semanage fcontext -a -t httpd_sys_rw_content_t '/apps/web/logs/nginx/nginx-error.log'
                                                        # restorecon -v '/apps/web/logs/nginx/nginx-error.log'


check status of the selinux and disable selinux temporary
```
sestatus
#disable selinux
setenforce 0
setenforce Permissive
```
These methods above will only work until the next reboot, therefore to disable SELinux permanently, move to the next section.

To permanently disable SELinux, use your favorite text editor to open the file /etc/sysconfig/selinux

By default, the SELinux configuration does not allow NGINX to connect to remote HTTP, FastCGI, or other servers, 
semanage  port -l | grep http_port_t

 semanage port -a -t http_port_t -p tcp 8000
ValueError: Port tcp/8000 already defined

modify


