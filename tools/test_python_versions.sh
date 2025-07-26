#! /bin/bash

while getopts 'l' flag; do
	case "${flag}" in
		l) local=true;;
	esac
done

source .venv310/bin/activate
echo "Testing in python 3.10"
if [ "$local" = true ]; then
	python spherapy/scripts/test_basic_functionality.py
else
	spherapy-test-basic
fi

source .venv311/bin/activate
echo "Testing in python 3.11"
if [ "$local" = true ]; then
	python spherapy/scripts/test_basic_functionality.py
else
	spherapy-test-basic
fi

source .venv312/bin/activate
echo "Testing in python 3.12"
if [ "$local" = true ]; then
	python spherapy/scripts/test_basic_functionality.py
else
	spherapy-test-basic
fi

deactivate