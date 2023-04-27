#/usr/bin/env python3

# 50 cores
# We distribute to 10 cores, on each one we run the main program with the size 10 times, obtaining 100 data points
# Then we get an average over all the results

from mpi4py import MPI
import subprocess
import numpy as np

from typing import Tuple

from datetime import datetime

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

repeats = 10
trieslimit = 100

basearr = np.linspace(start = 1, stop = 9, num = 9, endpoint = True, dtype = int)
arr = np.concatenate((basearr, basearr * 10, basearr * 100, basearr * 1000))
# arr = np.array([1,2,3,4,5,6,7,8,9,10])

print(arr)

numpts = len(arr)
timings_xgemm = np.zeros(shape = (repeats, numpts))
timings_naive = np.zeros(shape = (repeats, numpts))

for midx in range(numpts):
    m = arr[midx]

    print(rank, m)
    
    iteration = 0
    tries = 0
    while iteration < repeats:
        try:
            tries += 1
            if tries > trieslimit:
                timings_xgemm[iteration, midx] = np.nan
                timings_naive[iteration, midx] = np.nan
                iteration += 1
                continue

            if m <= 1000:
                timings_xgemm[iteration, midx], timings_naive[iteration, midx] = [int(x) for x in subprocess.check_output(["./main", f"{m}"]).strip().split(b"\n")]
            else:
                timings_xgemm[iteration, midx] = int(subprocess.check_output(["./main", f"{m}", "1"]).strip())
                timings_naive[iteration, midx] = np.nan
            iteration += 1
            
        except subprocess.CalledProcessError as e:
            print(rank, "Error: ", e)
            pass


recvBuffer_xgemm = None
recvBuffer_naive = None
if rank == 0:
    recvBuffer_xgemm = np.empty([size, repeats, numpts])
    recvBuffer_naive = np.empty([size, repeats, numpts])

comm.Gather(timings_xgemm, recvBuffer_xgemm, root = 0)
comm.Gather(timings_naive, recvBuffer_naive, root = 0)

if rank == 0: 
    def calculate(recvBuffer) -> Tuple[np.ndarray, np.ndarray]:
        collapsed = recvBuffer.reshape((repeats * size, numpts)).T
        avg    = np.nanmean(collapsed, axis = 1)
        stddev = np.nanstd(collapsed, axis = 1)

        return (avg, stddev)

    def patchNaN(num):
        if np.isnan(num):
            return ""
        return num

    xgemm = calculate(recvBuffer_xgemm)
    naive = calculate(recvBuffer_naive)

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d-%H:%M:%S")
    with open(f"slurm/scaling_result_{dt_string}.tsv", "w") as f:
        f.write("# m\tXGEMM\t\tNaive\n")
        f.write("# m\tRuntime(us)\tstddev(us)\tRuntime(us)\tstddev(us)\n")
        for midx in range(numpts):
            f.write(f"{arr[midx]}\t{xgemm[0][midx]}\t{xgemm[1][midx]}\t{patchNaN(naive[0][midx])}\t{patchNaN(naive[1][midx])}\n")