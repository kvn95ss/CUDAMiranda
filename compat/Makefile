# Standalone build for cudaMiranda (replaces the old CUDA SDK common.mk).
# Needs CUDA <= 11.x: the code uses legacy texture references, removed in CUDA 12.
NVCC  ?= /usr/bin/nvcc
ARCH  ?= -gencode arch=compute_86,code=compute_86
SRC   := ../src
FLAGS := -O3 -I. -I$(SRC) -I$(SRC)/H -w
NVFLAGS := $(FLAGS) $(ARCH) -Wno-deprecated-declarations

CFILES  := RNAlib_utils.c fold_vars.c energy_par.c fold.c params.c
CCFILES := output.cpp swat.cpp
CUFILES := miranda.cu scan.cu ExpString.cu pairs.cu seqio.cu thermo.cu utils.cu
OBJS := $(addprefix obj/,$(CFILES:.c=.c.o) $(CCFILES:.cpp=.cpp.o) $(CUFILES:.cu=.cu.o))

CUDAMiranda: $(OBJS)
	$(NVCC) $(ARCH) -o $@ $^ -lm
obj/%.c.o: $(SRC)/%.c | obj
	gcc $(FLAGS) -c $< -o $@
obj/%.cpp.o: $(SRC)/%.cpp | obj
	$(NVCC) $(NVFLAGS) -x cu -c $< -o $@
obj/%.cu.o: $(SRC)/%.cu | obj
	$(NVCC) $(NVFLAGS) -c $< -o $@
obj:
	mkdir -p obj
clean:
	rm -rf obj CUDAMiranda
