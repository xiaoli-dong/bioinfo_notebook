
# Install and configure nginx web server
yum module list nginx
yum install nginx

firewall-cmd --permanent --add-port={80/tcp,443/tcp}
firewall-cmd --reload
systemctl enable nginx
systemctl start nginx
systemctl status nginx
