# Problem 1 Implementing a matrix class

## 1.a
See C++ code and resulting `main`

## 1.b
`main` is designed to take 1 commandline argument, which is the size $m$ of the matrix to be copied, and outputs the time taken in microsecond to run that copy operation.

The python script `scaling.py` runs `main` with different $m$ values to obtain the time taken. For each value of $m$, `main m` was ran 80 times (20 iterations per CPU, 4 CPUs). The time is then aggregated over the 80 iterations to obtain the average and standard deviation.

This data is then plotted and fitted using `plot_scaling.gnu` to obtain the exponent $\alpha = 1.82$.

See runtime plot at `scaling.pdf`.