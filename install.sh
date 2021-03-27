#!/bin/sh

# Script for installation of PRISM-games on a clean install of Ubuntu

set -e # Abort if one of the commands fail
set -x # Print commands as they are executed

# Install dependencies: make/gcc/Java/git
sudo apt-get -y update
sudo apt -y install make gcc g++ default-jdk git

# Install Python (only needed for testing (prism-auto) currently)
sudo apt -y install python3 python3-pip
pip3 install matplotlib pandas scipy numpy

# Download the latest development version from GitHub
# git clone https://github.com/prismmodelchecker/prism-games.git

# Find the Java directory
export JAVA_HOME=`prism-games/prism/src/scripts/findjavac.sh | sed 's/\/bin\/javac//'`

# Download/build PPL + dependencies (needs Java 8 it seems)
# Comment out the section below if multi-obj support not needed
sudo apt -y install openjdk-8-jdk libgmp-dev m4
export JAVA8_HOME=`find /usr/lib/jvm -name 'java-8*'`

mkdir -p lib
cd lib
wget http://www.bugseng.com/products/ppl/download/ftp/releases/1.2/ppl-1.2.tar.gz
tar xfz ppl-1.2.tar.gz
cd ppl-1.2
./configure --enable-interfaces=java --with-java=$JAVA8_HOME
make
sudo make install
cd ..

# Download and build Yices 2.6
# Commented out for now since not enabled
# sudo apt -y install gperf
# wget https://yices.csl.sri.com/releases/2.6.0/yices-2.6.0-src.tar.gz
# tar xvfz yices-2.6.0-src.tar.gz
# cd yices-2.6.0
# ./configure
# make
# sudo make install
# cd ../../
# cp /usr/local/lib/libyices.so.2.6 prism-games/prism/lib

# Download and build Z3 4.8.7
# Commented out for now since z3 is included
# sudo apt -y install python
# wget https://github.com/Z3Prover/z3/archive/z3-4.8.7.tar.gz
# tar xvfz z3-4.8.7.tar.gz
# cd z3-z3-4.8.7
# python scripts/mk_make.py --java
# cd build
# make
# sudo make install
# cd ..
# cp /usr/lib/libz3.so prism-games/prism/lib

# Compile PRISM-games and run two tests
# (should ultimately display: "Testing result: PASS")
cd src/prism-games/prism
make PPL_DIR=/usr/local/lib
make test testppl

cd ../../prism-fruit/Games-DQL/
make PPL_DIR=/usr/local/lib


