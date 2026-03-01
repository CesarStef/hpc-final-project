#!/bin/bash

echo "Weak scaling: multinode scalability study with constant workload per resource"

N_STEPS=1000
TASKS_PER_NODE=8
OMP_NUM_THREADS=14
CPUS_PER_TASK=${OMP_NUM_THREADS}
LOCAL_SIZE=15000
GRID_SIZE_X=${LOCAL_SIZE}
GRID_SIZE_Y=${LOCAL_SIZE}

for NODES in 1 2 4 8 16; do
    TOTAL_PROCESS=$((NODES * TASKS_PER_NODE))

    GRID_SIZE_X=$((NODES * LOCAL_SIZE))

    if [ $NODES -eq 2 ]; then
        GRID_SIZE_X=$(( (LOCAL_SIZE) * 2 ))
        GRID_SIZE_Y=$(( LOCAL_SIZE ))
    elif [ $NODES -eq 4 ]; then
        GRID_SIZE_X=$(( LOCAL_SIZE * 2 ))
        GRID_SIZE_Y=$(( LOCAL_SIZE * 2 ))
    elif [ $NODES -eq 8 ]; then
        GRID_SIZE_X=$(( LOCAL_SIZE * 4 ))
        GRID_SIZE_Y=$(( LOCAL_SIZE * 2 ))
    elif [ $NODES -eq 16 ]; then
        GRID_SIZE_X=$(( LOCAL_SIZE * 4 ))
        GRID_SIZE_Y=$(( LOCAL_SIZE * 4 ))
    fi

    JOB_NAME="weak_scale_${NODES}n_${TOTAL_PROCESS}t"

    sbatch --nodes=${NODES} \
           --ntasks=${TOTAL_PROCESS} \
           --ntasks-per-node=${TASKS_PER_NODE} \
           --cpus-per-task=${CPUS_PER_TASK} \
           --job-name=${JOB_NAME} \
           --export=ALL,GRID_SIZE_X=${GRID_SIZE_X},GRID_SIZE_Y=${GRID_SIZE_Y},N_STEPS=${N_STEPS},OMP_NUM_THREADS=${OMP_NUM_THREADS},JOB_NAME=${JOB_NAME},TOTAL_PROCESS=${TOTAL_PROCESS} \
           run.sh

    echo "Submitting job with ${NODES} nodes, ${TOTAL_PROCESS} total tasks, grid size ${GRID_SIZE_X}x${GRID_SIZE_Y}"
done

echo "All Weak Scaling jobs submitted."