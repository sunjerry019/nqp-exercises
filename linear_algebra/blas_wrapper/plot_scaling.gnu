#!/usr/bin/env gnuplot

set term pdfcairo color size 6in, 5in
set output "scaling.pdf"

set title "Time taken to copy a square matrix of size $m$"
set xlabel "Size $m$"
set ylabel "Time Taken (microseconds)"

set mxtics
set mytics
set samples 10000

set grid mxtics xtics ytics mytics

set xrange [0:9200]

plot "slurm/scaling_result_20230425-18:05:20.tsv" u 1:2:3 with yerrorbars title "Data Points"