#!/bin/bash

# Run the 1D Linear modes problem for MHD. This script runs the 'slow', 
# 'alfven', and 'fast' waves for LINEAR and WENO5 reconstruction scheme at 4
# resolutions.

source ./utils.sh

# Make directory for mhdmodes1d outputs
MHDMODES1D_OUTPUT_DIR="$TESTS_OUTPUT_DIR/mhdmodes1d"
make_output_dir $MHDMODES1D_OUTPUT_DIR

# Path to problem directory
PROB="mhdmodes1d"
PROB_DIR="$HARM_DIR/prob/$PROB"

# Track current directory as we go about running the problem at different
# resolutions and then analyzing it
CURRENT_DIR=""

RECONS=("LINEAR" "WENO")
DIMENSIONS=(1 2)
MODES=("slow" "alfven" "fast")
RESOLUTIONS=(64 128 256 512)

for rec in "${RECONS[@]}"; do
    echo "Reconstruction scheme: $rec"
        
    RECON_DIR="$MHDMODES1D_OUTPUT_DIR/${rec}"
    make_output_dir $RECON_DIR

    for d in "${DIMENSIONS[@]}"; do
        echo "Problem along dimension: $d"
        
        DIM_DIR="$RECON_DIR/${d}dim"
        make_output_dir $DIM_DIR

        for m in "${MODES[@]}"; do
            echo "Mode: $m"

            MODES_DIR="$DIM_DIR/$m"
            make_output_dir $MODES_DIR

            for res in "${RESOLUTIONS[@]}"; do
                echo "Resolution: $res"

                RES_DIR="$MODES_DIR/$res"
                make_output_dir $RES_DIR

                # build harm
                build_harm $RES_DIR $PROB
                # copy over runtime par file
                cp $PROB_DIR/param.dat $RES_DIR
                # edit compile time params and build harm again
                if [[ "$d" == 1 ]]; then
                    set_res_x1 $RES_DIR $res
                    set_res_x2 $RES_DIR 1
                fi
                if [[ "$d" == 2 ]]; then
                    set_res_x1 $RES_DIR 1
                    set_res_x2 $RES_DIR $res
                fi
                set_recon_scheme $RES_DIR $rec
                build_harm $RES_DIR $PROB
                # edit runtime params
                # set_idim $d
                # set_mode $m
                # run harm
                # run_harm
            done
        done
    done
done