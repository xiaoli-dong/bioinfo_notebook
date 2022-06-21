



You first check out for the name of the package you want to remove:

dpkg --list

Then remove the given package

sudo apt-get remove package_name

Purge any related code

sudo apt-get purge package_name

Then Autoremove

sudo apt-get autoremove

Finally, do a clean so you check everything is correctly removed

sudo apt-get clean
