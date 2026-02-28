#!/bin/bash
#SBATCH --mem=0
#SBATCH --partition dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH -t 00:10:00
#SBATCH --exclusive

EXEC=./parallel

# =======================================================
module purge
module load openmpi/4.1.6--gcc--12.2.0

export OMP_PLACES=cores # Bind threads to cores, so that each thread gets its own core. This can improve performance by reducing contention for shared resources.
export OMP_PROC_BIND=close # Bind threads close to each other, which can improve cache locality and reduce communication overhead between threads.

if [[ -z "${TOTAL_PROCESS}" ]]; then
    export TOTAL_PROCESS=1
fi

if [[ -z "${OMP_NUM_THREADS}" ]]; then
    export OMP_NUM_THREADS=1
fi

if [[ ${TOTAL_PROCESS} -eq 1 ]]; then
    ${EXEC} -n ${N_STEPS} -x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -p 1
else
    mpirun -np ${TOTAL_PROCESS} ${EXEC} -n ${N_STEPS} -x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -p 1
fi