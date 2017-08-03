CFLAGS=-g3 -O0 -Wall -Werror -DUSE_FFT -fopenmp -pthread
OPTIM_CFLAGS=-g3 -O0 -Wall -Werror -DUSE_FFT -fopenmp -pthread
CXXFLAGS=-std=gnu++11 ${CFLAGS}
OPTIM_CXXFLAGS=-std=gnu++11 ${OPTIM_CFLAGS}

HEADERS=$(wildcard *.h) $(wildcard spqlios/*.h)

OBJS= \
spqlios/spqlios-fft-fma.o \
spqlios/spqlios-ifft-fma.o \
spqlios/spqlios-fft-impl.o \
spqlios/fft_processor_spqlios.o \
spqlios/lagrangehalfc_impl_fma.o \
test_parallel_multiplications.o \
poc_CircuitBootstrapping.o \
poc_karatsuba_optim.o

EXE=test_parallel_multiplications

%.o: %.s ${HEADERS}
	gcc -c ${OPTIM_CFLAGS} -o $@ $< 

%.o: %.cpp ${HEADERS}
	g++ -c ${CXXFLAGS} -o $@ $< 

%_optim.o: %.cpp ${HEADERS}
	g++ -c ${OPTIM_CXXFLAGS} -o $@ $< 

test_parallel_multiplications: ${OBJS}
	g++ ${CXXFLAGS} ${OBJS} -o ${EXE}

all: ${OBJS} ${EXE}

clean: 
	rm ${OBJS} ${EXE}
