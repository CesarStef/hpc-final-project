#!/bin/bash

NODES=1
N_TASKS_PER_NODE=1
TOTAL_PROCESS=1
N_STEPS=1000
GRID_SIZE_X=15000
GRID_SIZE_Y=15000

for OMP_THREADS in 1 2 4 8 16 32 56 84 112; do
    JOB_NAME="test_thread_${OMP_THREADS}_stefano_cattonar"
    export OMP_NUM_THREADS=${OMP_THREADS}

    sbatch --export=ALL,GRID_SIZE_X=${GRID_SIZE_X},GRID_SIZE_Y=${GRID_SIZE_Y},N_STEPS=${N_STEPS},JOB_NAME=${JOB_NAME},TOTAL_PROCESS=${TOTAL_PROCESS} --nodes=${NODES} --ntasks-per-node=${N_TASKS_PER_NODE} --cpus-per-task=${OMP_THREADS} --job-name=${JOB_NAME} run.sh
done