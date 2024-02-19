#!/bin/bash

DIR=$HOME/collision-melissa/msolve-app
TEST=$DIR/testapp/test.sh
BUMPER=$DIR/app/linux-x64/BumperTest

EXE=$BUMPER
#EXE=$TEST

$DIR/c-wrapper/build/heatc $EXE $1

