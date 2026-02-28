#!/bin/bash
gcc -Iinclude -Wall -march=native src/stencil_serial.c -o serial