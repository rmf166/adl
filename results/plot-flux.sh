#!/bin/sh
rm *.pdf
rm plot.p
echo 'set autoscale' >> plot.p
echo 'unset logscale' >> plot.p
echo 'unset label' >> plot.p
echo 'set size square' >> plot.p
echo 'set xtic auto' >> plot.p
echo 'set ytic auto' >> plot.p
echo 'set title "Problem 1"' >> plot.p
echo 'set xlabel "Position" enhanced' >> plot.p
echo 'set ylabel "Cell-Average Scalar Flux" enhanced' >> plot.p
echo 'set yr [0:40]' >> plot.p
echo 'set xr [0:5]' >> plot.p
echo 'plot "numave-p1-5-LD.dat"   using 1:2 title "LD"  with points pointtype 1 ps 1 lc rgb "red", \'  >> plot.p
echo '     "numave-p1-5-LC.dat"   using 1:2 title "LC"  with points pointtype 2 ps 1 lc rgb "blue", \' >> plot.p
echo '     "numave-p1-500-LD.dat" using 1:2 title "REF" with lines linetype        1 lc rgb "black"'   >> plot.p
echo 'set key left top' >> plot.p
echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
echo 'set output "plot-p1-ave.pdf"' >> plot.p
echo 'replot' >> plot.p
echo 'set terminal x11' >> plot.p
gnuplot plot.p
rm plot.p
echo 'set autoscale' >> plot.p
echo 'unset logscale' >> plot.p
echo 'unset label' >> plot.p
echo 'set size square' >> plot.p
echo 'set xtic auto' >> plot.p
echo 'set ytic auto' >> plot.p
echo 'set title "Problem 1"' >> plot.p
echo 'set xlabel "Position" enhanced' >> plot.p
echo 'set ylabel "Edge Scalar Flux" enhanced' >> plot.p
echo 'set yr [0:40]' >> plot.p
echo 'set xr [0:5]' >> plot.p
echo 'plot "numedg-p1-5-LD.dat"   using 1:2 title "LD"  with points pointtype 1 ps 1 lc rgb "red", \'  >> plot.p
echo '     "numedg-p1-5-LC.dat"   using 1:2 title "LC"  with points pointtype 2 ps 1 lc rgb "blue", \' >> plot.p
echo '     "numedg-p1-500-LD.dat" using 1:2 title "REF" with lines linetype        1 lc rgb "black"'   >> plot.p
echo 'set key left top' >> plot.p
echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
echo 'set output "plot-p1-edg.pdf"' >> plot.p
echo 'replot' >> plot.p
echo 'set terminal x11' >> plot.p
gnuplot plot.p
rm plot.p
echo 'set autoscale' >> plot.p
echo 'unset logscale' >> plot.p
echo 'unset label' >> plot.p
echo 'set xtic auto' >> plot.p
echo 'set ytic auto' >> plot.p
echo 'set title "Problem 2"' >> plot.p
echo 'set xlabel "Position" enhanced' >> plot.p
echo 'set ylabel "Cell-Average Scalar Flux" enhanced' >> plot.p
echo 'set yr [0:2]' >> plot.p
echo 'set xr [0:20]' >> plot.p
echo 'plot "numave-p2-20-LD.dat"   using 1:2 title "LD"  with points pointtype 1 ps 1 lc rgb "red", \'  >> plot.p
echo '     "numave-p2-20-LC.dat"   using 1:2 title "LC"  with points pointtype 2 ps 1 lc rgb "blue", \' >> plot.p
echo '     "numave-p2-2000-LD.dat" using 1:2 title "REF" with lines linetype        1 lc rgb "black"'   >> plot.p
echo 'set key left top' >> plot.p
echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
echo 'set output "plot-p2-ave.pdf"' >> plot.p
echo 'replot' >> plot.p
echo 'set terminal x11' >> plot.p
gnuplot plot.p
rm plot.p
echo 'set autoscale' >> plot.p
echo 'unset logscale' >> plot.p
echo 'unset label' >> plot.p
echo 'set xtic auto' >> plot.p
echo 'set ytic auto' >> plot.p
echo 'set title "Problem 2"' >> plot.p
echo 'set xlabel "Position" enhanced' >> plot.p
echo 'set ylabel "Edge Scalar Flux" enhanced' >> plot.p
echo 'set yr [0:2]' >> plot.p
echo 'set xr [0:20]' >> plot.p
echo 'plot "numedg-p2-5-LD.dat"   using 1:2 title "LD"  with points pointtype 1 ps 1 lc rgb "red", \'  >> plot.p
echo '     "numedg-p2-5-LC.dat"   using 1:2 title "LC"  with points pointtype 2 ps 1 lc rgb "blue", \' >> plot.p
echo '     "numedg-p2-500-LD.dat" using 1:2 title "REF" with lines linetype        1 lc rgb "black"'   >> plot.p
echo 'set key left top' >> plot.p
echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
echo 'set output "plot-p2-edg.pdf"' >> plot.p
echo 'replot' >> plot.p
echo 'set terminal x11' >> plot.p
gnuplot plot.p
pdfunite plot-p1-*.pdf plot-p2-*.pdf final.pdf
