#!/bin/sh

# remove installed library and dependencies

echo "Uninstalling library libwmgcr from /usr/lib"
if cd cpp; then
    if make uninstall; then
        echo "Library removed sucesfully"
    else
        echo >&2 "Error occured while removing library"
    fi
    echo -n "Removing compiled files from /cpp directory"
    if make clean; then
        echo "Compiled files removed sucesfully"
    fi
    echo -n "Remove armadillo library? (Y/n) "
    read answer
    # check if libarmadillo-dev is to be removed
    if [ "$answer" != "${answer#[Yy]}" ] ;then
        apt-get remove libarmadillo-dev
    fi
    # check if little cms is to be removed
    echo -n "Remove liblcms2-dev library? (Y/n) "
    read answer
    if [ "$answer" != "${answer#[Yy]}" ] ;then
        apt-get remove liblcms2-dev
    fi
else
    echo >&2 "Cpp directory is missing. Aborting"
    exit 1
fi