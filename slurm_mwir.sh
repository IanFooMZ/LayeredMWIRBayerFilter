#!/bin/bash

#SBATCH --time=144:00:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=8
#SBATCH --qos=normal
#SBATCH --mem-per-cpu=30G
#SBATCH --mail-user=ifoo@caltech.edu
#SBATCH --mail-type=END

source activate fdtd

xvfb-run --server-args="-screen 0 1280x1024x24" python LayeredMWIRBridgesBayerFilterOptimization.py 5 > stdout_mwir_a5.log 2> stderr_mwir_a5.log

exit $?
