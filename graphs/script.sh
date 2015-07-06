#!/bin/bash
g++ -c plot.cpp
g++ plot.o -o plot
./plot $1.plot	plotcmd	$1.pri	
gnuplot < plotcmd
latex plot.tex
dvips plot.dvi
gv plot.ps
ps2pdf14 plot.ps $1.pdf

