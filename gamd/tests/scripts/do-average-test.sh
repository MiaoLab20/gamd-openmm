#!/bin/bash

if ! command -v PyReweighting-1D.py &> /dev/null
then
    echo "You need to add the PyReweighting programs to your path prior to running this script."
    exit
fi

if [ "$#" -lt 3 ]; then
    echo "Usage do-average-test.sh location-of-cmd-directory run-type output-directory"
    exit 2
fi


GAMD_Directory=`realpath $1`
RUN_TYPE=$2
OUTPUT_BASE=$3
QUICK=$4
DEBUG=$5

mkdir $OUTPUT_BASE

./run-test.py $RUN_TYPE $OUTPUT_BASE/1/ $QUICK $DEBUG
echo ""
./run-test.py $RUN_TYPE $OUTPUT_BASE/2/ $QUICK $DEBUG
echo ""
./run-test.py $RUN_TYPE $OUTPUT_BASE/3/ $QUICK $DEBUG

COMPARISON_APP=`realpath ./tools/create-test-comparison-graphics.py`

cd  $OUTPUT_BASE/; $COMPARISON_APP $GAMD_Directory 1/ 2/ 3/; mv 1D-Phi.png $RUN_TYPE-1D-Phi.png; mv 1D-Psi.png $RUN_TYPE-1D-Psi.png


