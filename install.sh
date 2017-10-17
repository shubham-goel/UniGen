#!/bin/sh
cd ugen2/
autoreconf -f -i
mkdir build
cd build/
../configure
make -j2
sudo make install
