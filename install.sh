#!/bin/sh

# Check if the pkg-config available in the system. Neaded to find dependencies
if ! [ -x "$(command -v pkg-config)" ]; then
  echo >&2 'Error: pkg-config is not installed. Aborting.'
  exit 1
fi

echo "Check and install dependencies"

install_library(){
if pkg-config --exists $1; then
    echo "$1 library found!"
else 
    echo -n "$1 library missing! Install? (y/n)"
    read answer
    if [ "$answer" != "${answer#[Yy]}" ] ;then
        apt-get install lib$1-dev
    else
        echo >&2 "No dependency ($1) on the system. Aborting."
        exit 1
    fi
fi 
}

# Installing libraries if needed
install_library lapack
install_library blas
install_library armadillo
install_library lcms2

# compile and save library
echo "Compiling amd linking"
if cd cpp; then
    if make; then
        echo "Library compiled and linked sucesfully!"
    else
        echo >&2 "Error occured while compiling/linking. Aborting"
        exit 1
    fi
    echo -n "Install linked library libwmgcr.so to /usr/lib? (Y/n) "
    read answer
    if [ "$answer" != "${answer#[Yy]}" ] ;then
        make install
    else
        echo >&2 "Library is neaded for the on the system. Aborting."
        exit 1
    fi
else
    echo >&2 "Cpp directory is missing. Aborting"
    exit 1
fi