reset
set xlabel "car types"
set ylabel "time (s)"
set title "Time (4 households, maxage 25 years) with AD and Newton's method"
set term pdf
set output "plot1ad.pdf"
plot "data1ad.dat" u 2:($1==4?(stringcolumn(4) eq "C"?$5:1/0):1/0) title "Futhark c" with lines, \
     "" u 2:($1==4?(stringcolumn(4) eq "M"?$5:1/0):1/0) title "Futhark multicore" with lines, \
     "" u 2:($1==4?(stringcolumn(4) eq "O"?$5:1/0):1/0) title "Futhark opencl" with lines, \
     "matlabj.dat" u 2:($1==4?(stringcolumn(4) eq "BJ"?$5:1/0):1/0) title "Matlab (no AD)" with lines