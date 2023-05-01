#!/usr/bin/env gnuplot

# set term pdfcairo color size 6in, 5in
# set output "scaling.pdf"

set term pngcairo color size 12in, 10in
set output "scaling.png"

set title "Time taken to multiply square matrices of order $m$"
set xlabel "Size $m$"
set ylabel "Time Taken (microseconds)"

set mxtics
set mytics
set samples 10000

set grid mxtics xtics ytics mytics

f(x) = A*(x**B) + C
g(x) = D*(x**E) + F

fit f(x) "slurm/scaling_result_20230428-00:51:04.tsv" u 1:2:3 yerrors via A,B,C
fit g(x) "slurm/scaling_result_20230428-00:51:04.tsv" u 1:4:5 yerrors via D,E,F

set key top right spacing 1.5

set xrange [0:7200]
set yrange [0:1.2e8]

titelf = "XGEMM Fit: (".gprintf("%.5f", A).")m^{".gprintf("%.5f", B)."} + (".gprintf("%.5f", C).")"
titelg = "Naive Fit: (".gprintf("%.5f", D).")m^{".gprintf("%.5f", E)."} + (".gprintf("%.5f", F).")"
plot \
    "slurm/scaling_result_20230428-00:51:04.tsv" u 1:2:3 with yerrorbars title "XGEMM", \
    f(x) title titelf, \
    "slurm/scaling_result_20230428-00:51:04.tsv" u 1:4:5 with yerrorbars title "Naive", \
    g(x) title titelg

    # f(x) title titel, \
