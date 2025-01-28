#!/bin/bash

if [ "$1" == "" ]; then
    echo "Compiles and links the current project. This file should be located in the main project directory."
    echo ""
    echo "syntax:"
    echo "        m {build_type}"
    echo ""
    echo "input:"
    echo "       {build_type} - either optimised (o) or debug (d)" 
    exit
fi

# Pre-compile
rm -rf ${PROJECTDIR}/build
rm -rf ${PROJECTDIR}/bin
PROJECTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
let NPROCS=$(nproc | sed 's/[^0-9]//g')/2 # FLOOR DIVISION HERE

# Compile + Link
if [ "$1" == "o" ]; then
    cmake -B ${PROJECTDIR}/build/ -S ${PROJECTDIR}/. -DCMAKE_BUILD_TYPE=Release
elif [ "$1" == "d" ]; then
    cmake -B ${PROJECTDIR}/build/ -S ${PROJECTDIR}/. -DCMAKE_BUILD_TYPE=Debug
else
    echo "Invalid argument"
    exit
fi
cmake --build ${PROJECTDIR}/build/ -j $NPROCS
