#!/bin/bash

cd $HOME/collision-melissa/
if [ "$#" -eq 1 ]; then
				echo "Building files"
				./rebuild-cwrapper.sh
fi
echo "Deleting previous output."
find ./output -mindepth 1 \( -type d -name 'Data' -prune \) -o \( -exec rm -rf {} + \)
echo "Launching Melissa."
nohup melissa-launcher --config_name ./config/config_oar.json > run.out 2> run.err &
