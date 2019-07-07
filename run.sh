#!/bin/bash
clear
rm [1-9]* stuff.txt
mpicxx \-o main main.cpp && echo "Compilation Completed"
#mpiexec -np 4 ./main stuff.txt 13 77      | 78
mpiexec -np 4 ./main stuff.txt 34 1