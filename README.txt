QECONVERGE
=================================

Python code.

Before running DFT calculations, the ecutwfc and ecutrho values need to be chosen.  They are minimised such that the values are within a given converged limit.  They are also minimised because using an arbitrarily large cutoff would affect computer resources and computational time.

A suitable k-point mesh and smearing must also be determined.

Qeconverge takes an input file and a template pwscf file.  The position of the atoms can be varied slightly.  It starts at a specified ecutwfc and uses ecutrho = 4.0 x ecutwfc.  The value of ecutwfc (and ecutrho) are increased until the energy and forces converge.  The next stage is to reduce ecutrho while remaining in the convergence threshold (fixing ecutwfc).  Finally ecutrho is fixed and the code attemts to reduce ecutwfc further.

Further calculations are run and a 2D map ecutwfc on one axis and ecutrho on the other is mapped to visualise the convergence for the user.

Finally a 2D map with k-points on one axis and smearing on the other is created, to aid the user in choosing a k-point and smearing values.



Example 1:

examples/example1

The run file packages the python files into one called qeconverge.py - it then runs the python program.  System settings are stored in environment variables in the run.sh bash file:

export OMP_NUM_THREADS=1
export PROC_COUNT=4
export PWSCF_SCRATCH=/opt/scratch
export PWSCF_PP=/opt/pp
export PWSCF_CACHE=/opt/pwscf_cache
export PWSCF_BIN=/opt/qe/bin/pw.x


Example 2:

examples/example2

Same as example 1, but it runs the python file in its own directory.


Example 3:

examples/example3

This is similar to example 2, but the run.sh file is replaced with a batch.sh file - the names don't matter, but the batch.sh file contains more settings that are used on the BlueBEAR computer.  It is submitted to the job queue:

sbatch batch.sh

The settings may be modified to work on other supercomputers with different scheduling programs.  The python module and quantum espresso modules may beed to be loaded, although I use my own compiled quantum espresso.

batch.sh:

#!/bin/bash
#
#SBATCH --job-name=pwaleos
#SBATCH --output=jobout.txt
#SBATCH --account=readmsd02
#
#SBATCH --ntasks 20
#SBATCH --nodes 1
#SBATCH --time=360:00
#SBATCH --cpus-per-task=1
#
#SBATCH --get-user-env
#SBATCH --export=NONE
#
unset SLURM_EXPORT_ENV

module purge; module load bluebear
module load bear-apps/2018a
module load iomkl/2018a
module load Python/3.6.3-iomkl-2018a
module load matplotlib/2.1.1-iomkl-2018a-Python-3.6.3


# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
# Set the number of threads to 1
export OMP_NUM_THREADS=1
export PROC_COUNT=20
export PWSCF_SCRATCH=/scratch
export PWSCF_PP=/rds/homes/b/bxp912/pp
export PWSCF_CACHE=/rds/homes/b/bxp912/pwscf_cache
export PWSCF_BIN=/rds/homes/b/bxp912/apps/qe-6.3/bin/pw.x

python qeconverge.py input.in > result.txt





