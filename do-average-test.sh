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


GAMD_Directory=`realpath $1`
RUN_TYPE=$2
OUTPUT_BASE=$3

mkdir $OUTPUT_BASE

./run-test.py $RUN_TYPE $OUTPUT_BASE/1/
./run-test.py $RUN_TYPE $OUTPUT_BASE/2/
./run-test.py $RUN_TYPE $OUTPUT_BASE/3/

COMPARISON_APP=`realpath ./tools/create-test-comparison-graphics.py`

cd  $OUTPUT_BASE/; $COMPARISON_APP $GAMD_Directory 1/ 2/ 3/
