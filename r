#!/bin/bash

if [ "$1" == "" ]; then
    echo "Compiles and links the current project. This file should be located in the main project directory."
    echo ""
    echo "syntax:"
    echo "        r {np}"
    echo ""
    echo "input:"
    echo "       {np} - number of processors" 
    exit
fi

# Pre-execute
PROJECTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Execute
mpiexec -n $1 ${PROJECTDIR}/bin/HeatAdjoint
