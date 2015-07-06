#!/bin/bash
latex plot.tex
dvips plot.dvi
gv plot.ps
ps2pdf14 plot.ps $1.pdf

