#!/bin/bash
mpicc -Wall -Iinclude -fopenmp -O3 -march=native src/stencil_parallel.c -o parallel
