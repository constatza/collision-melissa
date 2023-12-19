#!/bin/sh
set -x

# Melissa will paste the `preprocessing_instructions`
export DOTNET_ROOT=$HOME/.dotnet
export PATH=$PATH:$DOTNET_ROOT:$DOTNET_ROOT/tools
../msolve-app/install-dotnet.sh

# melissa-launcher will find and replace 'melissa_set_env_file'
# automatically, do not change this line.


# User can set this part of the client script up automatically
# by ensuring that the keys in their `client_config` dictionary
# match the keywords below.  
# For example:
# melissa-launcher will search and replace "executable_command"
# automatically with the "executable_command" set in the 
# client_config file.

exec ../msolve-app/c-wrapper/build/heatc 5134 3 100 "$@"

