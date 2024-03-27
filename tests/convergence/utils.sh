#!/bin/bash

# A set of utility functions and paths common to all test problems.

# Path to iharm2d_v4 home directory
export HARM_DIR="../.." # local path because hard to specify universal global path
export HARM_DIR=$(realpath "$HARM_DIR") # now make it global
# Path to convergence tests output home directory
export TESTS_OUTPUT_DIR="$HARM_DIR/tests/convergence/tests_output"

# Check if TEST_OUTPUT_DIR exists; if not, make
if [ -d "$TESTS_OUTPUT_DIR" ]; then
    echo "Tests directory exists: $TESTS_OUTPUT_DIR"
else
    echo "Test directory does not exist. Creating: $TESTS_OUTPUT_DIR"
    mkdir -p "$TESTS_OUTPUT_DIR"
fi

# Function to create output directories if they don't exist
make_output_dir() {
    local OUTPUT_DIR="$1"
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Directory does not exist. Creating: $OUTPUT_DIR"
        mkdir -p "$OUTPUT_DIR"
    else
        echo "Directory exists. No action taken: $OUTPUT_DIR"
    fi
}

# Compile HARM
build_harm() {
    local SIM_DIR="$1"
    local PROB="$2"
    cd $SIM_DIR
    make -f $HARM_DIR/makefile PROB=$PROB
}

# Set X1 resolution
set_res_x1() {
    local SIM_DIR="$1"
    local N1="$2"
    sed -i -e "s/N1TOT [0-9]\+/N1TOT $N1/g" ${SIM_DIR}/build_archive/parameters.h
}

# Set X2 resolution
set_res_x2() {
    local SIM_DIR="$1"
    local N2="$2"
    sed -i -e "s/N2TOT [0-9]\+/N2TOT $N2/g" ${SIM_DIR}/build_archive/parameters.h

}

# Set reconstruction scheme
set_recon_scheme() {
    local SIM_DIR="$1"
    local RECON="$2"
    sed -i -e "s/RECONSTRUCTION [A-Z]\+/RECONSTRUCION $2/g" ${SIM_DIR}/build_archive/parameters.h
}

# Set mode dimension (for 1D linear modes)
set_idim() {
    local SIMDIR="$1"
    local idim="$2"
}