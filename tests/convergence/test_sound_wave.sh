#!/bin/bash

# Run the sound_wave test problem at 4 resolutions.

source ./utils.sh

# Make directory for mhdmodes1d outputs
SOUND_OUTPUT_DIR="$TESTS_OUTPUT_DIR/sound_wave"
make_output_dir $SOUND_OUTPUT_DIR

# Path to problem directory
PROB="sound_wave"
PROB_DIR="$HARM_DIR/prob/$PROB"

# Track current directory as we go about running the problem at different
# resolutions and then analyzing it
CURRENT_DIR=""

RECONS=("LINEAR" "WENO")
RESOLUTIONS=(64 128 256)

for rec in "${RECONS[@]}"; do
    echo "Reconstruction scheme: $rec"
    RECON_DIR="$SOUND_OUTPUT_DIR/${rec}"
    make_output_dir $RECON_DIR

    for res in "${RESOLUTIONS[@]}"; do
        echo "Resolution: $res"

        RES_DIR="$RECON_DIR/$res"
        make_output_dir $RES_DIR

        # build harm
        build_harm $RES_DIR $PROB
        # copy over runtime par file
        cp $PROB_DIR/param.dat $RES_DIR
        # edit compile time params and build harm again
        set_res_x1 $RES_DIR $res
        set_res_x2 $RES_DIR $res
        set_recon_scheme $RES_DIR $rec
        build_harm $RES_DIR $PROB
        # run harm
        run_harm $RES_DIR
    done

    # plot convergence
    cd $RECON_DIR
    cp $PROB_DIR/convergence_sound_wave.py ./
    IFS=,  # Set the Internal Field Separator to a comma
    res_comma_separated="${RESOLUTIONS[*]}"
    unset IFS
    python3 convergence_sound_wave.py --res=$res_comma_separated
done