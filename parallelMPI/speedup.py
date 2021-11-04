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


def generate_data(command:str, nprocesses:List[int], filename:str) -> List[float]:
    speedup = [1.0 for _ in range(len(nprocesses))]
    i = 0
    monothreadtime = 1.0
    with open(filename, "w") as f:
        for n in nprocesses:
            command_new = "mpirun -n {} {}".format(n, command)
            print("Lauching {}".format(command_new))
            process = subprocess.Popen(command_new.split(' '), stdout=subprocess.PIPE)
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
    subprocess.Popen("make").communicate()
    nprocesses = [i for i in range(1, 9)] # 9 is the max on my labtop
    nparticles = 100
    ntime = 5
    command1 = "./nbody_brute_force {} {}".format(nparticles, ntime)
    # command2 = "./nbody_barnes_hut {} {}".format(nparticles, ntime)

    subprocess.Popen("rm *.data", shell=True).communicate()
    speedup1 = generate_data(command1, nprocesses, "bruteforce.data")
    print(speedup1)
    # speedup2 = generate_data(command2, nprocesses, "barnes_hut.data")

    subprocess.Popen("rm *.png", shell=True).communicate()
    plot_speedup(nprocesses, {"MPI":speedup1}, "bruteforceMPIspeedup.png", "Brute force speedup")
    # plot_speedup(nprocesses, {"MPI":speedup2}, "barnes_hutMPIspeedup.png", "Barnes Hut speedup")

if __name__ == '__main__':
    main()