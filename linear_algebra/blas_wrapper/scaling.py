#/usr/bin/env python3

# 50 cores
# We distribute to 50 cores, on each one we run the main program with the size 20 times, obtaining 1000 data points
# Then we get an average over all the results
# Repeat until m = 5000

from mpi4py import MPI
import subprocess
import numpy as np

from datetime import datetime

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

start = 1
stop  = 10
repeats = 2

arr = np.linspace(start = start, stop = stop, num = (stop - start + 1), endpoint = True, dtype = int)
timings = np.zeros(shape = (len(arr), repeats))

shape = timings.shape
for midx in range(shape[0]):
    for iteration in range(shape[1]):
        m = arr[midx]
        timings[midx, iteration] = int(subprocess.check_output(["./main", f"{m}"]).strip())

recvBuffer = None
if rank == 0:
    recvBuffer = np.empty([size, shape[0], shape[1]])

comm.Gather(timings, recvBuffer, root = 0)

if rank == 0: 
    print(recvBuffer.shape)
    print(recvBuffer)

    avg = np.average(recvBuffer, axis = (0,2))

    # find stddev

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d-%H:%M:%S")
    with open(f"scaling_result_{dt_string}.tsv", "w") as f:
        f.write("# m\tRuntime\tstddev\n")
        for midx in range(shape[0]):
            f.write("")
