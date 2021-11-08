export PATH=/netfs/inf/trahay_f/mpich/bin/:$PATH

rm speedUpBruteforce.data
touch speedUpBruteforce.data

for nthreads in {1..12}
do
    for nproc in {1..8}
    do
        echo "compute with ${nthreads} threads and ${nproc} processes"
        duration=$(mpirun -np $nproc -machinefile machines ./nbody_brute_force 1000 1 $nthreads | grep "Simulation took" | cut -f3 -d ' ')
        echo "takes $duration s"
        echo "$nproc    $nthreads       $duration" >> speedUpBruteforce.data
    done
done
