#! /bin/python

from typing import Dict, List
import subprocess
import matplotlib.pyplot as plt


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
    i = 0
    with open(filename, "w") as f:
        for n in nthreads:
            print("Lauching {} with {} threads".format(command, n))
            process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE)
            output, _ = process.communicate()
            s = str(output.splitlines()[4])
            speedup[i] = float(s.split(' ')[2])
            print("time taken = {} s".format(speedup[i]))
            f.write("{}\t{}".format(n, speedup[i]))
            i += 1
    return speedup

def main():
    nthreads = [i for i in range(1, 16)]
    nparticles = 100
    ntime = 5
    command1 = "./nbody_brute_force {} {}".format(nparticles, ntime)
    command2 = "./nbody_barnes_hut {} {}".format(nparticles, ntime)

    speedup1 = generate_data(command1, nthreads, "bruteforce.data")
    speedup2 = generate_data(command2, nthreads, "barnes_hut.data")

    plot_speedup(nthreads, {"openMP":speedup1}, "bruteforceOMPspeedup.png", "Brute force speedup")
    plot_speedup(nthreads, {"openMP":speedup2}, "barnes_hutOMPspeedup.png", "Barnes Hut speedup")

if __name__ == '__main__':
    main()