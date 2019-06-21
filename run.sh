#!/bin/bash
clear
rm [1-9]*
mpicxx \-o main main.cpp
mpiexec -np 4 ./main stuff.txt 79 1
