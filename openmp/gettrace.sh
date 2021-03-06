#! /bin/bash
if [ $# -eq 0 ] ; then
	echo "usage: nbody_brute_force nparticules ntime (static or dynamic) nthreads"
    echo "usage: nbody_barnes_hut nparticules ntime (implem1 or implem2) nthreads"
	exit 0
fi
export OMP_NUM_THREADS=$5
make clean
make CC="eztrace_cc gcc" $1
eztrace -o /tmp -t omp ./$1 $2 $3 $4
eztrace_convert -o /tmp/eztrace_output_anquetil_trophime /tmp/${USER}_eztrace_log_rank_1
vite /tmp/eztrace_output_anquetil_trophime.trace