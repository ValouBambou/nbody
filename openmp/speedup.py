#! /bin/python

from typing import Dict, List
import subprocess
import matplotlib.pyplot as plt
import os


def plot_speedup(x_threads:List[int], y_dict:Dict[str, List[float]], filename:str, title:str) -> None:
    """plot and save fig of speedup given number of threads as x and labeled speedup as y"""
    plt.plot(x_threads, x_threads, label="linear")
    for key in y_dict:
        plt.plot(x_threads, y_dict[key], label=key)
    plt.title(title)
    plt.xlabel("number of threads")
    plt.ylabel("speedup")
    plt.legend()
    plt.savefig(filename)
    plt.show()


def generate_data(command:str, nthreads:List[int], filename:str) -> List[float]:
    speedup = [1.0 for _ in range(len(nthreads))]
    monothreadtime = 1.0
    i = 0
    with open(filename, "w") as f:
        for n in nthreads:
            print("Lauching {} with {} threads".format(command, n))
            process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, env={'OMP_NUM_THREADS':str(n), 'OMP_NESTED':'TRUE'})
            output, _ = process.communicate()
            s = str(output.splitlines()[4])
            if n == 1:
                monothreadtime = float(s.split(' ')[2])
            speedup[i] = monothreadtime / float(s.split(' ')[2])
            print("time taken = {} s".format(float(s.split(' ')[2])))
            f.write("{}\t{}\n".format(n, speedup[i]))
            i += 1
    return speedup


def main():
    nthreads = [i for i in range(1, 16)]
    nparticles = 2000
    ntime = 1

    subprocess.Popen("make").communicate()
    command_bfs = "./nbody_brute_force {} {} {}".format(nparticles, ntime, "static")
    command_bfd = "./nbody_brute_force {} {} {}".format(nparticles, ntime, "dynamic")
    command2 = "./nbody_barnes_hut {} {}".format(nparticles, ntime)

    subprocess.Popen("rm *.data", shell=True).communicate()
    speedup_bfs = generate_data(command_bfs, nthreads, "bruteforce_static.data")
    speedup_bfd = generate_data(command_bfd, nthreads, "bruteforce_dynamic.data")
    speedup2 = generate_data(command2, nthreads, "barnes_hut.data")

    subprocess.Popen("rm *.png", shell=True).communicate()
    plot_speedup(nthreads, {"OMPstatic":speedup_bfs, "OMPdynamic":speedup_bfd}, "bruteforceOMPspeedup.png", "Brute force speedup")
    plot_speedup(nthreads, {"openMP":speedup2}, "barnes_hutOMPspeedup.png", "Barnes Hut speedup")

if __name__ == '__main__':
    main()