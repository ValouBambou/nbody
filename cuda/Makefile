CC	= gcc

CUDA_HOME=/usr/local/cuda
CUDASDK=/usr/local/cuda/samples
CUDANVIDIA=/usr/lib64/nvidia

NVCC=$(CUDA_HOME)/bin/nvcc

NVCFLAGS= -I $(CUDA_HOME)/include -I$(CUDASDK)/common/inc

CFLAGS= -L $(CUDA_HOME)/lib64 -L $(CUDANVIDIA) -arch=sm_60 -lX11

VERBOSE	=
TARGET	= nbody_brute_force nbody_barnes_hut
OBJS	= ui.o xstuff.o nbody_tools.o nbody_alloc.o

#DISPLAY = -DDISPLAY
DISPLAY =

#DUMP = -DDUMP_RESULT
DUMP =

all: $(TARGET)

nbody_brute_force: nbody_brute_force.o $(OBJS)
	$(NVCC) $(VERBOSE) $(CFLAGS) -o $@ $< $(OBJS) 

nbody_barnes_hut: nbody_barnes_hut.o $(OBJS)
	$(NVCC) $(VERBOSE) $(CFLAGS) -o $@ $< $(OBJS)

nbody_brute_force.o: nbody_brute_force.cu
	$(NVCC) $(NVCFLAGS) $(CFLAGS) -c nbody_brute_force.cu $(VERBOSE) $(DISPLAY) $(DUMP)
%.o: %.c
	$(NVCC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY) $(DUMP)
clean:
	rm -f *.o $(TARGET)
