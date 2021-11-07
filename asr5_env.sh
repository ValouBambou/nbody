# MPI stuff
MPICH_ROOT=/usr/lib64//mpich/
export PATH=$PATH:$MPICH_ROOT/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MPICH_ROOT/lib
export MANPATH=/usr/share/man/mpich-x86_64/:${MANPATH}

#VampirTrace stuff
VT_ROOT=/opt/vampirtrace-5.14.4.1
export PATH=$PATH:$VT_ROOT/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VT_ROOT/lib

#PAPI stuff
#PAPI_ROOT=/opt/papi-5.5.1
#export PATH=$PATH:$PAPI_ROOT/bin
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PAPI_ROOT/lib

# CUDA stuff


export CUDA_HOME=/usr/local/cuda
export CUDA_INSTALL_PATH=${CUDA_HOME}
export CUDA_SDK_HOME=/usr/local/cuda/samples
export PATH=${CUDA_HOME}/bin:$PATH:.
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${CUDA_SDK_HOME}/common/lib/linux/x86>
export MANPATH=${CUDA_HOME}/man:$MANPATH

CUPTI_ROOT=/usr/local/cuda/extras/CUPTI/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUPTI_ROOT/lib64


#ViTE Stuff
VITE_ROOT=/opt/vite-1.2
export PATH=$PATH:$VITE_ROOT/bin

#EZTrace Stuff
EZTRACE_ROOT=/opt/eztrace-1.1-7
export PATH=$PATH:$EZTRACE_ROOT/bin
