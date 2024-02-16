#!/bin/bash
shopt -s extglob

cd $HOME/collision-melissa/
if [ "$#" -eq 1 ]; then
				echo "Building files"
				./rebuild-cwrapper.sh
fi

echo "Deleting previous output."
cd ./output
rm -rf !("InputFiles"|"SavedFiles")
cd ..

echo "Launching Melissa."
nohup melissa-launcher --config_name ./config/config_oar.json > run.out 2> run.err &
