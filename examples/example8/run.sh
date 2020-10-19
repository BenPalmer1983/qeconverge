#!/bin/bash

# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
# Set the number of threads to 1
export OMP_NUM_THREADS=1
export PROC_COUNT=4
export PWSCF_SCRATCH=/opt/scratch
export PWSCF_PP=/opt/pp
export PWSCF_CACHE=/opt/pwscf_cache
export PWSCF_BIN=/opt/qe/bin/pw.x


thisdir=$(pwd)
echo $thisdir
cd ../../
qeconverge=$(pwd)
echo $qeconverge
./package.sh
cd $thisdir
python3 $qeconverge/qeconverge.py input.in calc
#python3 $qeconverge/qeconverge.py input.in plot

