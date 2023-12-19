#!/bin/bash

cd $HOME/collision-melissa/
find ./output -mindepth 1 \( -type d -name 'Data' -prune \) -o \( -exec rm -rf {} + \)
melissa-launcher --config_name ./config/config_oar.json
