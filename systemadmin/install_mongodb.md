[Install MongoDB Community Edition on Red Hat](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-red-hat/)

open mongodb port for the remote connection

```
firewall-cmd --permanent --zone=public --add-port=27017/tcp
firewall-cmd --reload
firewall-cmd --list-ports
```
In windows, open mongodb compass and use the following connection string to connect to the remote server
```
mongodb://10.106.109.188:27017/
```
