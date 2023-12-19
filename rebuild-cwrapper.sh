#!/bin/bash

WORK_DIR=$HOME/collision-melissa/msolve-app/c-wrapper/
cd $WORK_DIR
rm -r ./build
mkdir build && cd build
cmake .. && make 

