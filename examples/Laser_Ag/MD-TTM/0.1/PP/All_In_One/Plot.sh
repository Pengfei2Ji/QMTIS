#! /bin/sh

python TTM_to_320.py

python Indiv_to_One.py

gnuplot < Dens.plt
convert Density.ps -rotate 90 Density.png

gnuplot < PressTens_x.plt
convert Presstens_x.ps -rotate 90 Presstens_x.png

gnuplot < Te.plt
convert Te.ps -rotate 90 Te.png

gnuplot < Tl.plt
convert Tl.ps -rotate 90 Tl.png
