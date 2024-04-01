#!/bin/bash

# Run Fishbone-Moncrief torus equilibrium test at 3 resolutions

source ./utils.sh

# Make directory for mhdmodes1d outputs
TORUS_OUTPUT_DIR="$TESTS_OUTPUT_DIR/torus_equilibrium"
make_output_dir $TORUS_OUTPUT_DIR

# Path to problem directory
PROB="torus_equilibrium"
PROB_DIR="$HARM_DIR/prob/$PROB"

# Track current directory as we go about running the problem at different
# resolutions and then analyzing it
CURRENT_DIR=""

RESOLUTIONS=(128 256 512)

for res in "${RESOLUTIONS[@]}"; do
    echo "Resolution: $res"

    RES_DIR="$TORUS_OUTPUT_DIR/$res"
    make_output_dir $RES_DIR

    # build harm
    build_harm $RES_DIR $PROB
    # copy over runtime par file
    cp $PROB_DIR/param.dat $RES_DIR
    # edit compile time params and build harm again
    set_res_x1 $RES_DIR $res
    set_res_x2 $RES_DIR $res
    build_harm $RES_DIR $PROB
    # run harm
    run_harm $RES_DIR
done
# plot convergence
cd $TORUS_OUTPUT_DIR
cp $PROB_DIR/convergence_torus.py ./
IFS=,  # Set the Internal Field Separator to a comma
res_comma_separated="${RESOLUTIONS[*]}"
unset IFS
python3 convergence_torus.py --res=$res_comma_separated