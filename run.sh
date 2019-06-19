#!/bin/bash
clear
rm [1-9]*
mpicxx \-o main main.cpp
mpiexec \-n 4 main stuff.txt 80 1
