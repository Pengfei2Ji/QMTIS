set terminal postscript color enhanced font 30 lw 3 size 18, 9

set xlabel "Time [ps]"
set ylabel "Position [65.3648 nm]"
set zlabel "Temperature [K]"

set output "Tl.ps"

set pm3d map
#set palette rgbformulae 22,13,-31
#set palette color
#set palette model XYZ rgbformulae 7,5,15
#set palette defined ( 0 "green", 1 "blue", 2 "red", 3 "orange" )
#set palette model XYZ rgbformulae 7,5,15

#set palette rgb 33, 13, 10

set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')


set xrange [20:160]
set yrange [0:14]

set xtics 20
set mxtics 5

set ytics 1
set mytics 5
set cbrange [0.0:2000]

#set pm3d map
#set contour surface
#set cntrparam levels discrete 0, 300, 1235, 2435, 8000
#set samples 2435
#set isosamples 2435


#set contour base


#set contour surface
#set cntrparam levels discrete 300, 1234.93, 2435, 5000
#set cntrlabel  format '%8.3g' font '7' start 5 interval 1000
#set cntrparam levels auto 1000
#set style data lines

splot "Gnuplot_Output.dat" using ($1*0.5):($2*13.85/320.00):($6) with p pointtype 5 pointsize 1 palette linewidth 100 notitle

#splot "Gnuplot_Den_Output.dat" using 1:2:3 with points palette pointsize 3 pointtype 7
#splot "Gnuplot_Den_Output.dat" using 1:2:3 with l palette pointsize 3 pointtype 7
