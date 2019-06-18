#!/bin/bash
clear
rm [1-9]*
mpicxx \-o main main.cpp
mpiexec \-n 1 main stuff.txt 66 1
