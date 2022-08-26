


# Uninstalled Ubuntu and debian packages with apt-get & dpkg
Get a List of Installed Packages:
```
dpkg --list
```

Remove package using apt-get
```
# remove the given package
sudo apt-get remove package_name

#Purge any related code
sudo apt-get purge package_name

# Then Autoremove
sudo apt-get autoremove

#Finally, do a clean so you check everything is correctly removed
sudo apt-get clean
```
