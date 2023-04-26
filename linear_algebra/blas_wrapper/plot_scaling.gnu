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

A = 2
B = 2
f(x) = A*(x**B) + C

fit f(x) "slurm/scaling_result_20230425-18:05:20.tsv" u 1:2:3 yerrors via A,B,C
# fit f(x) "slurm/scaling_result_20230425-18:05:20.tsv" u 1:2 via A

set key bottom right spacing 1.5

titel = "Fit: (".gprintf("%.5f", A).")m^{".gprintf("%.5f", B)."} + (".gprintf("%.5f", C).")"
plot f(x) title titel, \
    "slurm/scaling_result_20230425-18:05:20.tsv" u 1:2:3 with yerrorbars title "Data Points"