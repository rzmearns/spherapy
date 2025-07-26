#!/bin/bash

while getopts 'l' flag; do
	case "${flag}" in
		l) local=true;;
	esac
done

# Remove any old virtual environments
rm -rf .venv310
rm -rf .venv311
rm -rf .venv312

# create virtual env
virtualenv -p=`which python3.10` .venv310
virtualenv -p=`which python3.11` .venv311
virtualenv -p=`which python3.12` .venv312


# install
source .venv310/bin/activate
if [ "$local" = true ]; then
	pip install --editable .[dev]
else
	pip install git+ssh://git@gitlab.unimelb.edu.au/msl/libraries/spherapy.git
fi

source .venv311/bin/activate
if [ "$local" = true ]; then
	pip install --editable .[dev]
else
	pip install git+ssh://git@gitlab.unimelb.edu.au/msl/libraries/spherapy.git
fi

source .venv312/bin/activate
if [ "$local" = true ]; then
	pip install --editable .[dev]
else
	pip install git+ssh://git@gitlab.unimelb.edu.au/msl/libraries/spherapy.git
fi

deactivate