FLAGS=-DNDEBUG


all:
	g++ -O0 -g3 -march=native -Wall -std=gnu++11 code.cpp -o code -lntl

optim: 
	g++ ${FLAGS} -Ofast -march=native -fno-stack-protector -g3 -std=gnu++11 code.cpp -o code -lntl

bench_fft_gmp: bench_fft_gmp.cpp
	g++ -O0 -g3 -Wall -std=gnu++11 $< -o $@ -lmpfr -lgmp -lm

bench_fft_gmp_optim: bench_fft_gmp.cpp
	g++ ${FLAGS} -Ofast -march=native -fno-stack-protector -g3 -std=gnu++11 $< -o $@ -lmpfr -lgmp -lm


bench_fft_mpfr: bench_fft_mpfr.cpp
	g++ -O0 -g3 -Wall -std=gnu++11 $< -o $@ -lmpfr -lgmp -lm

bench_fft_mpfr_optim: bench_fft_mpfr.cpp
	g++ ${FLAGS} -Ofast -march=native -fno-stack-protector -g3 -std=gnu++11 $< -o $@ -lmpfr -lgmp -lm
