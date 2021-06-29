#!/bin/bash

if ! command -v PyReweighting-1D.py &> /dev/null
then
    echo "You need to add the PyReweighting programs to your path prior to running this script."
    exit
fi

if [ "$#" -ne 3 ]; then
    echo "Usage do-average-test.sh location-of-gamd-directory run-type output-directory"
    exit 2
fi

GAMD_Directory=$1
RUN_TYPE=$2
OUTPUT_BASE=$3

./run-test.py $RUN_TYPE $OUTPUT_BASE-1
./run-test.py $RUN_TYPE $OUTPUT_BASE-2
./run-test.py $RUN_TYPE $OUTPUT_BASE-3

./tools/create-test-comparison-graphics.py $GAMD_Directory $OUTPUT_BASE-1 $OUTPUT_BASE-2 $OUTPUT_BASE-3
