reset
set xlabel "car types"
set ylabel "time (s)"
set title "Time (4 households, maxage 25 years) with Newton's method"
set term pdf
set output "plot1j.pdf"
plot "data1j.dat" u 2:($1==4?(stringcolumn(4) eq "C"?$5:1/0):1/0) title "Futhark c" with lines, \
     "" u 2:($1==4?(stringcolumn(4) eq "M"?$5:1/0):1/0) title "Futhark multicore" with lines, \
     "" u 2:($1==4?(stringcolumn(4) eq "O"?$5:1/0):1/0) title "Futhark opencl" with lines, \
     "matlabj.dat" u 2:($1==4?(stringcolumn(4) eq "BJ"?$5:1/0):1/0) title "Matlab" with lines