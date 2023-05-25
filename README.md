pXar
====

pixel Xpert analysis &amp; readout

## Prerequisites

- [libusb-1.0]() 
    - can be installed with `apt install libusb-1.0-0-dev` on Ubuntu 
- [FTDI chip drivers](https://ftdichip.com/drivers/d2xx-drivers/) for Linux.
```shell
            # get the latest driver for you system and extract it
            curl -L -O "https://ftdichip.com/wp-content/uploads/2022/07/libftd2xx-x86_64-1.4.27.tgz"
            # The unpacker puts everything in the directory called release
            tar xvf libftd2xx-x86_64-1.4.27.tgz 
            # copy the libraries to where linker can find it and create sim links
            sudo cp release/build/lib* /usr/local/lib
            sudo ln -sf /usr/local/lib/libftd2xx.so.1.4.27 /usr/local/lib/libftd2xx.so
            # Copy the header files as well to the include dir where it can be found by your compiler
            sudo cp WinTypes.h ftd2xx.h /usr/local/include/
```

## Installation
```shell
git clone git@github.com:diamondIPP/pxar.git
cd pxar/
mkdir build
cd build/
cmake .. 
# optionally use ccmake or other cmake GUI to fix the problems
make install

