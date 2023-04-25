#/usr/bin/env python3

# 50 cores
# We distribute to 4 cores, on each one we run the main program with the size 20 times, obtaining 800 data points
# Then we get an average over all the results

from mpi4py import MPI
import subprocess
import numpy as np

from datetime import datetime

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

repeats = 20

basearr = np.linspace(start = 1, stop = 9, num = 9, endpoint = True, dtype = int)
arr = np.concatenate((basearr, basearr * 10, basearr * 100, basearr * 1000))

print(arr)

numpts = len(arr)
timings = np.zeros(shape = (repeats, numpts))

for midx in range(numpts):
    m = arr[midx]

    print(rank, m)
    
    for iteration in range(repeats):
        timings[iteration, midx] = int(subprocess.check_output(["./main", f"{m}"]).strip())

recvBuffer = None
if rank == 0:
    recvBuffer = np.empty([size, repeats, numpts])

comm.Gather(timings, recvBuffer, root = 0)

if rank == 0: 
    collapsed = recvBuffer.reshape((repeats * size, numpts)).T

    avg    = np.average(collapsed, axis = 1)
    stddev = np.std(collapsed, axis = 1)

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d-%H:%M:%S")
    with open(f"slurm/scaling_result_{dt_string}.tsv", "w") as f:
        f.write("# m\tRuntime(us)\tstddev(us)\n")
        for midx in range(numpts):
            f.write(f"{arr[midx]}\t{avg[midx]}\t{stddev[midx]}\n")