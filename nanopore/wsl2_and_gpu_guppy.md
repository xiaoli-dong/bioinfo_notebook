# How to setup WSL2 on windows 10, and enable gpu based guppy basecalling
Here is a list of references I used during my setup process

- [CUDA on WSL User Guide](https://docs.nvidia.com/cuda/wsl-user-guide/index.html)  
- [Nanopore Guppy GPU basecalling on Windows using WSL2](https://hackmd.io/@Miles/rkYKDHPsO)  
- [Nanopore software download](https://community.nanoporetech.com/downloads)  
- [How to mount windows network drives in wsl](https://www.public-health.uiowa.edu/it/support/kb48568/)

In order to enable the terminal tabs and panel split, you will need to install Windows Terminal through [Microsoft Store](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701?hl=en-au&gl=AU) or you can download it from [github](https://github.com/microsoft/terminal). Here is a reference page on [how to split, resize, switch between panels](https://docs.microsoft.com/en-us/windows/terminal/panes) after the installation

### Troubleshooting while I am doing setup

When I updated the system, I got the following error message

```
xdong@M691822:~$ sudo apt-get update
Get:1 http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64  InRelease [1575 B]
Hit:2 http://security.ubuntu.com/ubuntu focal-security InRelease
Err:1 http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64  InRelease
  The following signatures couldn't be verified because the public key is not available: NO_PUBKEY A4B469963BF863CC
Hit:3 http://archive.ubuntu.com/ubuntu focal InRelease
Hit:4 http://archive.ubuntu.com/ubuntu focal-updates InRelease
Hit:5 http://archive.ubuntu.com/ubuntu focal-backports InRelease
Reading package lists... Done
W: GPG error: http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64  InRelease: The following signatures couldn't be verified because the public key is not available: NO_PUBKEY A4B469963BF863CC
E: The repository 'http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64  InRelease' is not signed.
N: Updating from such a repository can't be done securely, and is therefore disabled by default.
N: See apt-secure(8) manpage for repository creation and user configuration details.
```

Here is [the reference](https://chrisjean.com/fix-apt-get-update-the-following-signatures-couldnt-be-verified-because-the-public-key-is-not-available/) I used to fix the problem.   
```
# add missing keys, run the following commands:
xdong@M691822:~$ sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys A4B469963BF863CC
```

