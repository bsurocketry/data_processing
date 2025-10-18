set title "corrected thrust"
set xlabel "seconds"
set ylabel "newtons"
plot "prior_data.txt" using 1:2 with lines title "prior",   \
     "post_data.txt" using 1:2 with lines title "corrected"
