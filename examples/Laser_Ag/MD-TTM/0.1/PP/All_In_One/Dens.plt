set terminal postscript color enhanced font 30 lw 3 size 18, 9

set xlabel "Time [ps]"
set ylabel "Position [65.3648 nm]"
set zlabel "Density [g/cm^3]"

set output "Density.ps"

set pm3d map
#set palette rgbformulae 22,13,-31
set palette color
#set palette model XYZ rgbformulae 7,5,15
#set palette defined ( 0 "green", 1 "blue", 2 "red", 3 "orange" )
#set palette model XYZ rgbformulae 7,5,15

#set palette rgb 33, 13, 10

#set palette defined ( 0 '#000090',\
#                      1 '#000fff',\
#                      2 '#0090ff',\
#                      3 '#0fffee',\
#                      4 '#90ff70',\
#                      5 '#ffee00',\
#                      6 '#ff7000',\
#                      7 '#ee0000',\
#                      8 '#7f0000')


set xrange [20:160]
set yrange [0:14]

set xtics 20
set mxtics 5

set ytics 1
set mytics 5

set cbrange [0.1:1.2]

#set ticslevel 
#set cntrparam levels 10
#set auto
#set zrange [-1.0:1.0]
#set style data lines
#set title "Density Contour"
#set parametric
#set view map

#set view map
#set size ratio .9

#set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back
#set object 1 rect fc rgb "black" fillstyle solid 1.0

splot "Gnuplot_Output.dat" using ($1*0.5):($2*13.85/320.00):($3) with p pointtype 5 pointsize 1 palette linewidth 100 notitle

#splot "Gnuplot_Den_Output.dat" using 1:2:3 with points palette pointsize 3 pointtype 7
#splot "Gnuplot_Den_Output.dat" using 1:2:3 with l palette pointsize 3 pointtype 7
