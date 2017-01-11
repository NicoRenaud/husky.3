#!/bin/bash

# input file
IN=./hsk.in

# output file
OUT=./

#EXEC
EXE=../husky

# compute the TE
$EXE $IN $OUT

# plot the TE
cd $OUT
gnuplot plot_results.dem

