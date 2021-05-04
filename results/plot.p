set autoscale
unset logscale
unset label
set xtic auto
set ytic auto
set title "Problem 2"
set xlabel "Position" enhanced
set ylabel "Edge Scalar Flux" enhanced
set yr [0:2]
set xr [0:20]
plot "numedg-p2-5-LD.dat"   using 1:2 title "LD"  with points pointtype 1 ps 1 lc rgb "red", \
     "numedg-p2-5-LC.dat"   using 1:2 title "LC"  with points pointtype 2 ps 1 lc rgb "blue", \
     "numedg-p2-500-LD.dat" using 1:2 title "REF" with lines linetype        1 lc rgb "black"
set key left top
set terminal pdfcairo enhanced color dashed
set output "plot-p2-edg.pdf"
replot
set terminal x11
