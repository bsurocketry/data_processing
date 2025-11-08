set title "values"
set xlabel "seconds"
set ylabel "newtons"
set y2label "pascals"
set y2tics
set y2range [0 to 4e6]
plot "prior_data.txt" using 1:2 with lines title "prior",   \
     "post_data.txt" using 1:2 with lines title "corrected", \
     "chamber_pressure.txt" using 1:2 with lines title "chamber pressure" axis x1y2
