#!/bin/bash

wd=$(pwd)
cd ../
./package.sh
cd $wd
cp ../pwscf_converge.py pwscf_converge.py
./run.sh
