set xlabel "x"
set ylabel "y"
set zlabel "z"

set term qt
splot 'data.txt' title "Data" lw 6, 'out.txt' title "Interpolation"
set term png
set output "example.png"
rep
set term qt
