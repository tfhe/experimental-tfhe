#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>

using namespace std;

typedef int32_t Torus32; //avant uint32_t
//typedef int64_t Torus64; //avant uint64_t



default_random_engine generator;
uniform_int_distribution<Torus32> uniformTorus32_distrib(INT32_MIN, INT32_MAX);
uniform_int_distribution<int> uniformInt_distrib(INT_MIN, INT_MAX);

/** sets the seed of the random number generator to the given values */
EXPORT void tfhe_random_generator_setSeed(uint32_t* values, int size) {
    seed_seq seeds(values, values+size);
    generator.seed(seeds);
}

// Gaussian sample centered in message, with standard deviation sigma
EXPORT Torus32 gaussian32(Torus32 message, double sigma){
    //Attention: all the implementation will use the stdev instead of the gaussian fourier param
    normal_distribution<double> distribution(0.,sigma); //TODO: can we create a global distrib of param 1 and multiply by sigma?
    double err = distribution(generator);
    return message + dtot32(err);
}



// from double to Torus32
EXPORT Torus32 dtot32(double d) {
    return int32_t(int64_t((d - int64_t(d))*_two32));
}
// from Torus32 to double
EXPORT double t32tod(Torus32 x) {
    return double(x)/_two32_double;
}


EXPORT Torus32 approxPhase(Torus32 phase, int Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t half_interval = interv/2; // begin of the first intervall
    uint64_t phase64 = (uint64_t(phase)<<32) + half_interval;
    //floor to the nearest multiples of interv
    phase64 -= phase64%interv;
    //rescale to torus32
    return int32_t(phase64>>32); 
}

EXPORT int modSwitchFromTorus32(Torus32 phase, int Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t half_interval = interv/2; // begin of the first intervall
    uint64_t phase64 = (uint64_t(phase)<<32) + half_interval;
    //floor to the nearest multiples of interv
    return phase64/interv;
}

EXPORT Torus32 modSwitchToTorus32(int mu, int Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t phase64 = mu*interv;
    //floor to the nearest multiples of interv
    return phase64>>32;
}































//allocate memory space for a LagrangeHalfCPolynomial
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial() {
    return (LagrangeHalfCPolynomial*) malloc(sizeof(LagrangeHalfCPolynomial));
}
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial_array(int nbelts) {
    return (LagrangeHalfCPolynomial*) malloc(nbelts*sizeof(LagrangeHalfCPolynomial));
}

//free memory space for a LweKey
EXPORT void free_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* ptr) {
    free(ptr);
}
EXPORT void free_LagrangeHalfCPolynomial_array(int nbelts, LagrangeHalfCPolynomial* ptr) {
    free(ptr);
}

//allocates and initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ new)
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial(const int N) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial();
    init_LagrangeHalfCPolynomial(obj,N);
    return obj;
}
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int nbelts, const int N) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial_array(nbelts);
    init_LagrangeHalfCPolynomial_array(nbelts,obj,N);
    return obj;
}

//destroys and frees the LagrangeHalfCPolynomial structure
//(equivalent of the C++ delete)
EXPORT void delete_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj) {
    destroy_LagrangeHalfCPolynomial(obj);
    free_LagrangeHalfCPolynomial(obj);
}
EXPORT void delete_LagrangeHalfCPolynomial_array(int nbelts, LagrangeHalfCPolynomial* obj) {
    destroy_LagrangeHalfCPolynomial_array(nbelts,obj);
    free_LagrangeHalfCPolynomial_array(nbelts,obj);
}

/** multiplication via direct FFT (it must know the implem of LagrangeHalfCPolynomial because of the tmp+1 notation */
EXPORT void torusPolynomialMultFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(result, tmp+2);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
    TorusPolynomial* tmpr = new_TorusPolynomial(N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(tmpr, tmp+2);
    torusPolynomialAddTo(result, tmpr);
    delete_TorusPolynomial(tmpr);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
EXPORT void torusPolynomialSubMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
    TorusPolynomial* tmpr = new_TorusPolynomial(N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(tmpr, tmp+2);
    torusPolynomialSubTo(result, tmpr);
    delete_TorusPolynomial(tmpr);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}

























// TorusPolynomial = 0
EXPORT void torusPolynomialClear(TorusPolynomial* result) {
    const int N = result->N;

    for (int i = 0; i < N; ++i) result->coefsT[i] = 0;
}

// TorusPolynomial = random
EXPORT void torusPolynomialUniform(TorusPolynomial* result) {
    const int N = result->N;
    Torus32* x = result->coefsT;

    for (int i = 0; i < N; ++i) 
    x[i] = uniformTorus32_distrib(generator);
}

// TorusPolynomial = TorusPolynomial
EXPORT void torusPolynomialCopy(
    TorusPolynomial* result, 
    const TorusPolynomial* sample) {
    assert(result!=sample);
    const int N = result->N;
    const Torus32* __restrict s = sample->coefsT;
    Torus32* __restrict r = result->coefsT;

    for (int i = 0; i < N; ++i) r[i] = s[i];
}

// TorusPolynomial + TorusPolynomial
EXPORT void torusPolynomialAdd(TorusPolynomial* result, const TorusPolynomial* poly1, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    assert(result!=poly1); //if it fails here, please use addTo
    assert(result!=poly2); //if it fails here, please use addTo
    Torus32* __restrict r = result->coefsT;
    const Torus32* __restrict a = poly1->coefsT;
    const Torus32* __restrict b = poly2->coefsT;

    for (int i = 0; i < N; ++i) 
    r[i] = a[i] + b[i];
}

// TorusPolynomial += TorusPolynomial
EXPORT void torusPolynomialAddTo(TorusPolynomial* result, const TorusPolynomial* poly2) {
    const int N = poly2->N;
    Torus32* r = result->coefsT;
    const Torus32* b = poly2->coefsT;

    for (int i = 0; i < N; ++i) 
    r[i] += b[i];
}


// TorusPolynomial - TorusPolynomial
EXPORT void torusPolynomialSub(TorusPolynomial* result, const TorusPolynomial* poly1, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    assert(result!=poly1); //if it fails here, please use subTo
    assert(result!=poly2); //if it fails here, please use subTo
    Torus32* __restrict r = result->coefsT;
    const Torus32* a = poly1->coefsT;
    const Torus32* b = poly2->coefsT;

    for (int i = 0; i < N; ++i) 
    r[i] = a[i] - b[i];
}

// TorusPolynomial -= TorusPolynomial
EXPORT void torusPolynomialSubTo(TorusPolynomial* result, const TorusPolynomial* poly2) {
    const int N = poly2->N;
    Torus32* r = result->coefsT;
    const Torus32* b = poly2->coefsT;

    for (int i = 0; i < N; ++i) 
    r[i] -= b[i];
}

// TorusPolynomial + p*TorusPolynomial
EXPORT void torusPolynomialAddMulZ(TorusPolynomial* result, const TorusPolynomial* poly1, int p, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    Torus32* r = result->coefsT;
    const Torus32* a = poly1->coefsT;
    const Torus32* b = poly2->coefsT;

    for (int i = 0; i < N; ++i) 
    r[i] = a[i] + p*b[i];
}

// TorusPolynomial += p*TorusPolynomial
EXPORT void torusPolynomialAddMulZTo(TorusPolynomial* result, const int p, const TorusPolynomial* poly2) {
    const int N = poly2->N;
    Torus32* r = result->coefsT;
    const Torus32* b = poly2->coefsT;

    for (int i = 0; i < N; ++i) r[i] += p*b[i];
}

// TorusPolynomial - p*TorusPolynomial
EXPORT void torusPolynomialSubMulZ(TorusPolynomial* result, const TorusPolynomial* poly1, const int p, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    Torus32* r = result->coefsT;
    const Torus32* a = poly1->coefsT;
    const Torus32* b = poly2->coefsT;

    for (int i = 0; i < N; ++i) r[i] = a[i] - p*b[i];
}

//result= (X^{a}-1)*source
EXPORT void torusPolynomialMulByXaiMinusOne(TorusPolynomial* result, int a, const TorusPolynomial* source){
    const int N=source->N;
    Torus32* out=result->coefsT;
    Torus32* in =source->coefsT; 

    assert(a>=0 && a<2*N);

    if (a<N) {
    for (int i=0;i<a;i++)//sur que i-a<0
        out[i]= -in[i-a+N]-in[i];
    for (int i=a;i<N;i++)//sur que N>i-a>=0
        out[i]= in[i-a]-in[i];
    } else {
    const int aa=a-N;
    for (int i=0;i<aa;i++)//sur que i-a<0
        out[i]= in[i-aa+N]-in[i];
    for (int i=aa;i<N;i++)//sur que N>i-a>=0
        out[i]= -in[i-aa]-in[i];
    }
}


//result= X^{a}*source
EXPORT void torusPolynomialMulByXai(TorusPolynomial* result, int a, const TorusPolynomial* source){
    const int N=source->N;
    Torus32* out=result->coefsT;
    Torus32* in =source->coefsT; 

    assert(a>=0 && a<2*N);
    assert(result != source);

    if (a<N) {
    for (int i=0;i<a;i++)//sur que i-a<0
        out[i]= -in[i-a+N];
    for (int i=a;i<N;i++)//sur que N>i-a>=0
        out[i]= in[i-a];
    } else {
    const int aa=a-N;
    for (int i=0;i<aa;i++)//sur que i-a<0
        out[i]= in[i-aa+N];
    for (int i=aa;i<N;i++)//sur que N>i-a>=0
        out[i]= -in[i-aa];
    }
}


// TorusPolynomial -= p*TorusPolynomial
EXPORT void torusPolynomialSubMulZTo(TorusPolynomial* result, int p, const TorusPolynomial* poly2) {
    const int N = poly2->N;
    Torus32* r = result->coefsT;
    const Torus32* b = poly2->coefsT;

    for (int i = 0; i < N; ++i) r[i] -= p*b[i];
}


// Norme Euclidienne d'un IntPolynomial
EXPORT double intPolynomialNormSq2(const IntPolynomial* poly){
    const int N = poly->N;
    int temp1 = 0;

    for (int i = 0; i < N; ++i)
    {
        int temp0 = poly->coefs[i]*poly->coefs[i];
        temp1 += temp0;
    }
    return temp1;
}

// Sets to zero
EXPORT void intPolynomialClear(IntPolynomial* poly){
    const int N = poly->N;
    for (int i = 0; i < N; ++i)
        poly->coefs[i]=0;
}

// Sets to zero
EXPORT void intPolynomialCopy(IntPolynomial* result, const IntPolynomial* source){
    const int N = source->N;
    for (int i = 0; i < N; ++i)
        result->coefs[i]=source->coefs[i];
}

/** accum += source */
EXPORT void intPolynomialAddTo(IntPolynomial* accum, const IntPolynomial* source) {
    const int N = source->N;
    for (int i = 0; i < N; ++i)
        accum->coefs[i]+=source->coefs[i];
}

/**  result = (X^ai-1) * source */
EXPORT void intPolynomialMulByXaiMinusOne(IntPolynomial* result, int ai, const IntPolynomial* source) {
    const int N=source->N;
    int* out=result->coefs;
    int* in =source->coefs; 

    assert(ai>=0 && ai<2*N);

    if (ai<N) {
    for (int i=0;i<ai;i++)//sur que i-a<0
        out[i]= -in[i-ai+N]-in[i];
    for (int i=ai;i<N;i++)//sur que N>i-a>=0
        out[i]= in[i-ai]-in[i];
    } else {
    const int aa=ai-N;
    for (int i=0;i<aa;i++)//sur que i-a<0
        out[i]= in[i-aa+N]-in[i];
    for (int i=aa;i<N;i++)//sur que N>i-a>=0
        out[i]= -in[i-aa]-in[i];
    }
}



// Norme infini de la distance entre deux TorusPolynomial
EXPORT double torusPolynomialNormInftyDist(const TorusPolynomial* poly1, const TorusPolynomial* poly2) {
    const int N = poly1->N;
    double norm = 0;

    // Max between the coefficients of abs(poly1-poly2)
    for (int i = 0; i < N; ++i){
        double r = abs(t32tod(poly1->coefsT[i] - poly2->coefsT[i]));
        if (r>norm) {norm = r;}
    }
    return norm;
}






// Norme 2 d'un IntPolynomial
EXPORT double intPolynomialNorm2sq(const IntPolynomial* poly) {
    const int N = poly->N;
    double norm = 0;

    for (int i = 0; i < N; ++i){
        double r = poly->coefs[i];
        norm += r*r;
    }
    return norm;
}

// Norme infini de la distance entre deux IntPolynomial
EXPORT double intPolynomialNormInftyDist(const IntPolynomial* poly1, const IntPolynomial* poly2) {
    const int N = poly1->N;
    double norm = 0;


    // Max between the coefficients of abs(poly1-poly2)
    for (int i = 0; i < N; ++i){
        double r = abs(poly1->coefs[i] - poly2->coefs[i]);
        if (r>norm) {norm = r;}
    }
    return norm;
}


