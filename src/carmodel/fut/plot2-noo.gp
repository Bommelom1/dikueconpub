reset
set xlabel "ntypes (households)"
set ylabel "time (s)"
set title "Time (35 car types, maxage 25 years)"
set term pdf
set output "plot2.pdf"
plot "data2.dat" u 1:($3==25?(stringcolumn(4) eq "C"?$5:1/0):1/0) title "Futhark c" with lines, \
     "" u 1:($3==25?(stringcolumn(4) eq "M"?$5:1/0):1/0) title "Futhark multicore" with lines, \
     "matlab2.dat" u 1:($3==25?(stringcolumn(4) eq "B"?$5:1/0):1/0) title "Matlab" with lines