#ifndef LAGRANGEHALFC_IMPL_SPQLIOS_H
#define LAGRANGEHALFC_IMPL_SPQLIOS_H

#include <cassert>
#include <cmath>
#include <cstdint>

class FFT_Processor_Spqlios {
    public:
    const int _2N;
    const int N;    
    const int Ns2;
    private:
    double* real_inout_direct;
    double* imag_inout_direct;
    double* real_inout_rev;
    double* imag_inout_rev;
    void* tables_direct;
    void* tables_reverse;
    
    public:
    FFT_Processor_Spqlios(const int N);

    void execute_reverse_int(double* res, const int* a);
    void execute_reverse_torus32(double* res, const int32_t* a);
    void execute_direct_torus32(int32_t* res, const double* a);
    void execute_reverse_torus64(double* res, const int64_t* a);
    void execute_direct_torus64(int64_t* res, const double* a);

    ~FFT_Processor_Spqlios();
};

extern FFT_Processor_Spqlios fftp1024;
extern FFT_Processor_Spqlios fftp2048;

extern "C" void LagrangeHalfCPolynomialAddMulASM(double* res, double* a, double* b, long Ns2);

#endif // LAGRANGEHALFC_IMPL_SPQLIOS_H
