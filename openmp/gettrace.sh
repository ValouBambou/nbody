#! /bin/bash

make clean
make CC="eztrace_cc gcc" $1
eztrace -o /tmp -t omp ./$1 $2 $3
eztrace_convert -o /tmp/eztrace_output /tmp/${USER}_eztrace_log_rank_1
vite /tmp/eztrace_output.trace