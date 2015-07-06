#!/bin/bash
rm format.tex
rm matrixexpr.txt
rm /home/jacobc/Desktop/Desktop1026/clustalw/format2.txt
 ./seq2exp -bs fas.txt -bi topbot2Int.txt -s rhoseq2.txt -p synmyout6 -e rhoexp2.tab -m factordts.wtmx -f factorexpdts2.tab -na 1 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -nrand 5 -et 7 -dc coreOPT0dc1.txt -du coreOPT0du1.txt -sa NEE.txt -e2 NEEexp.tab

