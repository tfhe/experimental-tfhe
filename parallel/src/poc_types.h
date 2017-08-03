#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>

#include "generic_utils.h"

#define k 1

//basic types
typedef int32_t Torus32; // x represents the Torus value x/2^32 mod 1
typedef int64_t Torus64; // x represents the Torus value x/2^64 mod 1

static const int64_t _two32 = INT64_C(1) << 32; // 2^32
inline Torus32 t64tot32(Torus64 x) {
    return int32_t(x/_two32);
}
inline Torus64 t32tot64(Torus32 x) {
    return int64_t(x)*_two32;
}




/*
struct PolynomialParameter32 {
    const int N; 
    void* const FFT_PREPROC;
};

struct PolynomialParameter64 {
    const int N;
    void* const FFT_PREPROC;
};
*/


struct Torus32Polynomial {
    Torus32* const coefs;
    // constructor 
    Torus32Polynomial(int N): coefs(new Torus32[N]) {}
    // destructor 
    ~Torus32Polynomial() { delete[] coefs; }
};
/*
// alloc and initialize
Torus32Polynomial* new_Torus32Polynomial(int N) { 
    Torus32Polynomial* obj = (Torus32Polynomial*) malloc(sizeof(Torus32Polynomial));
    new(obj) Torus32Polynomial(N);
    return obj;
}
Torus32Polynomial* new_Torus32Polynomial_array(int nbelts, int N) {
    Torus32Polynomial* obj = (Torus32Polynomial*) malloc(nbelts*sizeof(Torus32Polynomial));
    for (int i = 0; i < nbelts; i++) new(obj+i) Torus32Polynomial(N);
    return obj;
}
// destroy and free 
void delete_Torus32Polynomial(Torus32Polynomial* obj) { 
    obj->~Torus32Polynomial();
    free(obj);
}
void delete_Torus32Polynomial_array(int nbelts, Torus32Polynomial* obj) {
    for (int i = 0; i < nbelts; i++) (obj+i)->~Torus32Polynomial();
    free(obj);
}
*/




struct Torus64Polynomial {
    Torus64* const coefs;
    // constructor 
    Torus64Polynomial(int N): coefs(new Torus64[N]) {}
    // destructor 
    ~Torus64Polynomial() { delete[] coefs; }
};





struct IntPolynomial {
    int* const coefs;
    // constructor 
    IntPolynomial(int N): coefs(new int[N]) {}
    // destructor 
    ~IntPolynomial() { delete[] coefs; }
};



#ifdef USE_FFT
struct LagrangeHalfCPolynomial {
    double* const values;
    LagrangeHalfCPolynomial(int N): values(new double[N]) {}
    ~LagrangeHalfCPolynomial() {
        delete[] values;
    }
};
#else
struct LagrangeHalfCPolynomial {
    IntPolynomial* intPoly; 
    Torus64Polynomial* torus64Poly; 
    LagrangeHalfCPolynomial(int N) {
        intPoly=0;
        torus64Poly=0;
    }
    void clear() {
        if (intPoly) { delete intPoly; intPoly=0; }
        if (torus64Poly) { delete torus64Poly; torus64Poly=0;}
    }
    void setIntPoly(const IntPolynomial* a, const int N) {
        clear();
        intPoly = new IntPolynomial(N);
        for (int i=0; i<N; i++) intPoly->coefs[i]=a->coefs[i];
    }
    void setTorus64Poly(const Torus64Polynomial* a, const int N) {
        clear();
        torus64Poly = new Torus64Polynomial(N);
        for (int i=0; i<N; i++) torus64Poly->coefs[i]=a->coefs[i];
    }
    void setZeroTorus64Poly(const int N) {
        clear();
        torus64Poly = new Torus64Polynomial(N);
        for (int i=0; i<N; i++) torus64Poly->coefs[i]=0;
    }
    ~LagrangeHalfCPolynomial() {
    }
};
#endif



struct LweSample32 {
    Torus32* const a;
    Torus32* const b; //b = &a[n]
    // constructor 
    LweSample32(int n): a(new Torus32[n+1]), b(&a[n]) {}
    // destructor 
    ~LweSample32() { delete[] a; }    
};






struct LweSample64 {
    Torus64* const a;
    Torus64* const b; //b = &a[n]
    // constructor 
    LweSample64(int n): a(new Torus64[n+1]), b(&a[n]) {}
    // destructor 
    ~LweSample64() { delete[] a; }  
};





struct TLweSample32 {
    Torus32Polynomial* const a;
    Torus32Polynomial* const b; //b = &a[k]
    // constructor 
    TLweSample32(int N): a(new_array1<Torus32Polynomial>(k+1,N)), b(&a[k]) {}
    // destructor 
    ~TLweSample32() { delete_array1<Torus32Polynomial>(a); }  
};





struct TLweSample64 {
    Torus64Polynomial* const a;
    Torus64Polynomial* const b; //b = &a[k]
    // constructor 
    TLweSample64(int N): a(new_array1<Torus64Polynomial>(k+1,N)), b(&a[k]) {}
    // destructor 
    ~TLweSample64() { delete_array1<Torus64Polynomial>(a); }  
};





struct TLweSampleFFT {
    LagrangeHalfCPolynomial* const a;
    LagrangeHalfCPolynomial* const b; //b = &a[k]
    // constructor 
    TLweSampleFFT(int N): a(new_array1<LagrangeHalfCPolynomial>(k+1,N)), b(&a[k]) {}
    // destructor 
    ~TLweSampleFFT() { delete_array1<LagrangeHalfCPolynomial>(a); }  
};








struct TGswSample32 {
    TLweSample32** const samples; //samples[k+1][ell]
    TLweSample32* const allsamples;
    // constructor 
    TGswSample32(int l, int N): 
        samples(new_array2<TLweSample32>((k+1),l,N)),
        allsamples(samples[0]) {}
    // destructor 
    ~TGswSample32() { 
        delete_array2<TLweSample32>(samples);
    }      
};





struct TGswSample64 {
    TLweSample64** const samples; //samples[k+1][ell]
    TLweSample64* const allsamples;
    // constructor 
    TGswSample64(int l, int N): 
        samples(new_array2<TLweSample64>((k+1),l,N)),
        allsamples(samples[0]) {}
    // destructor 
    ~TGswSample64() { 
        delete_array2<TLweSample64>(samples);
    }     
};





struct TGswSampleFFT {
    TLweSampleFFT** const samples; //samples[k+1][ell]
    TLweSampleFFT* const allsamples;
    // constructor 
    TGswSampleFFT(int l, int N): 
        samples(new_array2<TLweSampleFFT>((k+1),l,N)),
        allsamples(samples[0]) {}
    // destructor 
    ~TGswSampleFFT() { 
        delete_array2<TLweSampleFFT>(samples);
    }  
};















class Globals {
    public:
        //primary parameters
        static const int n_lvl0;
        static const int n_lvl1;
        static const int n_lvl2;
        static const int bgbit_lvl1;
        static const int ell_lvl1;
        static const int bgbit_lvl2;
        static const int ell_lvl2;
        static const double bkstdev_lvl2;
        static const double ksstdev_lvl10;
        static const int kslength_lvl10;
        static const int ksbasebit_lvl10;
        static const double ksstdev_lvl21;
        static const int kslength_lvl21;
        static const int ksbasebit_lvl21;

        //deduced parameters
        int t_lvl0; // t_lvl0 = kslength_lvl10 * ksbasebit_lvl10 ?
        int t_lvl1; // t_lvl1 = kslength_lvl21 * ksbasebit_lvl21 
        int N_lvl1; // N_lvl1 = n_lvl1
        int N_lvl2; // N_lvl2 = n_lvl2
        // Offset
        uint64_t torusDecompOffset;
        // Buffer
        uint64_t* torusDecompBuf; 
        
        //secret keys (and their polynomial interpretation)
        int* key_lvl0;
        int* key_lvl1;
        IntPolynomial* Key_lvl1;
        int* key_lvl2;
        IntPolynomial* Key_lvl2;
        
        //cloud keys
        //PreKS
        LweSample32*** preKS; // preKS[n_lvl1][kslength_lvl10][ksbase_lvl10]
        //circuit bootstrapping keys
        TGswSample64* bk; // bk[n_lvl0]
        TGswSampleFFT* bkFFT; // bkFFT[n_lvl0]
        //privKS
        TLweSample32**** privKS;

        Globals();
}; 


