#! /bin/python

from typing import Dict, List
import subprocess
import matplotlib.pyplot as plt


def plot_speedup(ncores: Dict[str, List[int]], y_dict: Dict[str, List[float]], filename: str, title: str) -> None:
    """plot and save fig of speedup given number of cores as x and labeled speedup as y"""
    plt.plot(ncores, ncores, label="linear")
    for key in y_dict:
        plt.plot(ncores[key], y_dict[key], label=key)
    plt.title(title)
    plt.xlabel("number of cores")
    plt.ylabel("speedup")
    plt.legend()
    plt.savefig(filename)
    plt.show()


def main():
    max_nthreads = 12
    max_nproc = 8
    data_matrix = [[0.0 for _ in range(max_nthreads)]
                   for _ in range(max_nproc)]
    with open("speedUpBruteforce.data") as f:
        for line in f.readlines():
            s = line.split('    ')
            i = int(s[0]) - 1
            j = int(s[1]) - 1
            value = float(s[2])
            data_matrix[i][j] = value

    d = {"{} MPI processes".format(i+1): [data_matrix[i][j] / data_matrix[0][0]
                                          for j in range(max_nthreads)] for i in range(max_nproc)}
    ncores = {"{} MPI processes".format(
        i): [i*j for j in range(1, max_nthreads+1)] for i in range(1, max_nproc)}
    plot_speedup(ncores, d, "MPI_OMP_speedup.png", "Speedup fonction de nthreads * nproccess")


if __name__ == '__main__':
    main()
