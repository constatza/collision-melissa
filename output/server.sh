#!/bin/sh
set -x

# Melissa will paste the `preprocessing_instructions`
. /home/catzarakis/spack/share/spack/setup-env.sh
spack load melissa-api
spack load py-melissa-core

# the remainder of this file should be left untouched. 
# melissa-launcher will find and replace values in 
# curly brackets (e.g) with 
# the proper values.


echo "DATE                      =$(date)"
echo "Hostname                  =$(hostname -s)"
echo "Working directory         =$(pwd)"
echo ""
echo $PYTHONPATH

set -e

exec melissa-server --project_dir /home/catzarakis/collision-melissa/./config --config_name config_oar
