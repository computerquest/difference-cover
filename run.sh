#!/bin/bash
clear
rm [1-9]* stuff.txt
mpicxx \-o main main.cpp
#mpiexec -np 4 ./main stuff.txt 13 77
mpiexec -np 2 ./main stuff.txt 65 1