#! /bin/sh

tar -xvf vcflib.tar.gz
cd vcflib/tabixpp/htslib
make
cd ../..
make
cd ..

tar -xvf sdsl-lite.tar.gz
cd sdsl-lite
./install.sh "$(pwd)"
cd ..
