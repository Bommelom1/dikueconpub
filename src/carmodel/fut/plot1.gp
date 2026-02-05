reset
set xlabel "car types"
set ylabel "time (s)"
set title "Time (4 households, maxage 20 years)"
set term pdf
set output "plot1.pdf"
plot "data1.dat" u 2:($1==4?(stringcolumn(4) eq "C"?$5:1/0):1/0) title "Futhark c" with lines, \
     "" u 2:($1==4?(stringcolumn(4) eq "M"?$5:1/0):1/0) title "Futhark multicore" with lines, \
     "" u 2:($1==4?(stringcolumn(4) eq "O"?$5:1/0):1/0) title "Futhark opencl" with lines, \
     "matlab.dat" u 2:($1==4?(stringcolumn(4) eq "B"?$5:1/0):1/0) title "Matlab" with lines