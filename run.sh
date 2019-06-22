#!/bin/bash
clear
rm [1-9]*
mpicxx \-o main main.cpp
#mpiexec -np 4 ./main stuff.txt 13 77
mpiexec -np 4 ./main stuff.txt 22 1