#!/bin/sh

# Script for installation of PRISM-games on a clean install of Ubuntu

set -e # Abort if one of the commands fail
set -x # Print commands as they are executed

# Install dependencies: make/gcc/Java/git
sudo apt-get -y update
sudo apt -y install make gcc g++ openjdk-8-jdk git

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
cd ../../

# Compile PRISM-games and run two tests
# (should ultimately display: "Testing result: PASS")
pwd
cd src/prism-games/prism

make PPL_DIR=/usr/local/lib
make test testppl

cd ../../prism-fruit/Games-DQL/
make PPL_DIR=/usr/local/lib


