#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>


using namespace std;

typedef int32_t Torus32; //avant uint32_t
//typedef int64_t Torus64; //avant uint64_t

 
TLweParams::TLweParams(int N, int k, double alpha_min, double alpha_max):
    N(N),
    k(k),
    alpha_min(alpha_min),
    alpha_max(alpha_max),
    extracted_lweparams(N*k,alpha_min, alpha_max){}
TLweParams::~TLweParams() {}

TLweKey::TLweKey(const TLweParams* params):params(params){
    key = new_IntPolynomial_array(params->k,params->N);   
}
TLweKey::~TLweKey() {delete_IntPolynomial_array(params->k,key);}

TLweSample::TLweSample(const TLweParams* params): k(params->k) {
    a = new_TorusPolynomial_array(k+1, params->N);
    b = a+k;
    current_variance = 0;
}
TLweSample::~TLweSample() {delete_TorusPolynomial_array(k+1, a);}

TLweSampleFFT::TLweSampleFFT(const TLweParams* params, LagrangeHalfCPolynomial* arr, double current_variance): k(params->k) {
    a = arr;
    b = a+k;
    current_variance = 0;
}
TLweSampleFFT::~TLweSampleFFT() {}









/* TLWE functions */
EXPORT void tLweKeyGen(TLweKey* result){
    const int N = result->params->N;
    const int k = result->params->k;
    uniform_int_distribution<int> distribution(0,1);

    for (int i = 0; i < k; ++i)
        for (int j = 0; j < N; ++j)
            result->key[i].coefs[j] = distribution(generator);
}

EXPORT void tLweSymEncryptZero(TLweSample* result, double alpha, const TLweKey* key){
    const int N = key->params->N;
    const int k = key->params->k;
    
    for (int j = 0; j < N; ++j)
        result->b->coefsT[j] = gaussian32(0, alpha);   
    
    for (int i = 0; i < k; ++i) {
    torusPolynomialUniform(&result->a[i]);
    torusPolynomialAddMulR(result->b, &key->key[i], &result->a[i]);
    }

    result->current_variance = alpha*alpha;
}

EXPORT void tLweSymEncrypt(TLweSample* result, TorusPolynomial* message, double alpha, const TLweKey* key){
    const int N = key->params->N;
    
    tLweSymEncryptZero(result, alpha, key);

    for (int j = 0; j < N; ++j)
        result->b->coefsT[j] += message->coefsT[j];   
}

EXPORT void tLweSymEncryptT(TLweSample* result, Torus32 message, double alpha, const TLweKey* key){
    tLweSymEncryptZero(result, alpha, key);

    result->b->coefsT[0] += message;
}



EXPORT void tLwePhase(TorusPolynomial* phase, const TLweSample* sample, const TLweKey* key){
    const int k = key->params->k;

    torusPolynomialCopy(phase, sample->b); // phi = b

    for (int i = 0; i < k; ++i)
        torusPolynomialSubMulR(phase, &key->key[i], &sample->a[i]);
}


EXPORT void tLweApproxPhase(TorusPolynomial* message, const TorusPolynomial* phase, int Msize, int N){
    for (int i = 0; i < N; ++i) message->coefsT[i] = approxPhase(phase->coefsT[i], Msize);
}


EXPORT void tLweSymDecrypt(TorusPolynomial* result, const TLweSample* sample, const TLweKey* key, int Msize){
    tLwePhase(result, sample, key);
    tLweApproxPhase(result, result, Msize, key->params->N);
}


EXPORT Torus32 tLweSymDecryptT(const TLweSample* sample, const TLweKey* key, int Msize){
    TorusPolynomial* phase = new_TorusPolynomial(key->params->N);

    tLwePhase(phase, sample, key);
    Torus32 result = approxPhase(phase->coefsT[0], Msize);

    delete_TorusPolynomial(phase);
    return result;
}


EXPORT void tLweClear(TLweSample* result, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i < k; ++i) torusPolynomialClear(&result->a[i]);
    torusPolynomialClear(result->b);
    result->current_variance = 0.;
}


EXPORT void tLweCopy(TLweSample* result, const TLweSample* sample, const TLweParams* params){
    const int k = params->k;
    const int N = params->N;

    for (int i = 0; i <= k; ++i) 
        for (int j = 0; j < N; ++j)
            result->a[i].coefsT[j] = sample->a[i].coefsT[j];
    
    result->current_variance = sample->current_variance;
}



EXPORT void tLweNoiselessTrivial(TLweSample* result, const TorusPolynomial* mu, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i < k; ++i) torusPolynomialClear(&result->a[i]);
    torusPolynomialCopy(result->b, mu);
    result->current_variance = 0.;
}

EXPORT void tLweNoiselessTrivialT(TLweSample* result, const Torus32 mu, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i < k; ++i) torusPolynomialClear(&result->a[i]);
    torusPolynomialClear(result->b);
    result->b->coefsT[0]=mu;
    result->current_variance = 0.;
}

EXPORT void tLweAddTo(TLweSample* result, const TLweSample* sample, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i < k; ++i) 
    torusPolynomialAddTo(&result->a[i], &sample->a[i]);
    torusPolynomialAddTo(result->b, sample->b);
    result->current_variance += sample->current_variance;
}

EXPORT void tLweSubTo(TLweSample* result, const TLweSample* sample, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i < k; ++i) 
    torusPolynomialSubTo(&result->a[i], &sample->a[i]);
    torusPolynomialSubTo(result->b, sample->b);
    result->current_variance += sample->current_variance; 
}

EXPORT void tLweAddMulTo(TLweSample* result, int p, const TLweSample* sample, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i < k; ++i) 
    torusPolynomialAddMulZTo(&result->a[i], p, &sample->a[i]);
    torusPolynomialAddMulZTo(result->b, p, sample->b);
    result->current_variance += (p*p)*sample->current_variance;
}

EXPORT void tLweSubMulTo(TLweSample* result, int p, const TLweSample* sample, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i < k; ++i) 
    torusPolynomialSubMulZTo(&result->a[i], p, &sample->a[i]);
    torusPolynomialSubMulZTo(result->b, p, sample->b);
    result->current_variance += (p*p)*sample->current_variance; 
}

EXPORT void tLweAddMulRTo(TLweSample* result, const IntPolynomial* p, const TLweSample* sample, const TLweParams* params){
    const int k = params->k;
    
    for (int i = 0; i <= k; ++i)
       torusPolynomialAddMulR(result->a+i, p, sample->a+i); 
    result->current_variance += intPolynomialNormSq2(p)*sample->current_variance; 
}


//mult externe de X^ai-1 par bki
EXPORT void tLweMulByXaiMinusOne(TLweSample* result, int ai, const TLweSample* bk, const TLweParams* params){
    const int k=params->k;
    for(int i=0;i<=k;i++)
        torusPolynomialMulByXaiMinusOne(&result->a[i],ai,&bk->a[i]);
}



EXPORT void init_TLweKey(TLweKey* obj, const TLweParams* params) {
    new(obj) TLweKey(params);
}
EXPORT void destroy_TLweKey(TLweKey* obj) {
    (obj)->~TLweKey();
}

EXPORT void init_TLweSample(TLweSample* obj, const TLweParams* params) {
    new(obj) TLweSample(params);
}
EXPORT void destroy_TLweSample(TLweSample* obj) {
    (obj)->~TLweSample();
}








































/* Pas sure qu'on ait besoin de tout cela!!! */

EXPORT void init_TLweSampleFFT(TLweSampleFFT* obj, const TLweParams* params) {
    //a is a table of k+1 polynomials, b is an alias for &a[k]
    const int k = params->k;
    LagrangeHalfCPolynomial* a = new_LagrangeHalfCPolynomial_array(k+1, params->N);
    double current_variance = 0;
    new(obj) TLweSampleFFT(params, a, current_variance);
}


EXPORT void destroy_TLweSampleFFT(TLweSampleFFT* obj) {
    const int k = obj->k;
    delete_LagrangeHalfCPolynomial_array(k+1, obj->a);
    obj->~TLweSampleFFT();
}



EXPORT void tLweToFFTConvert(TLweSampleFFT* result, const TLweSample* source, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i <= k; ++i)
    TorusPolynomial_ifft(result->a+i,source->a+i);
    result->current_variance=source->current_variance;
}



EXPORT void tLweFromFFTConvert(TLweSample* result, const TLweSampleFFT* source, const TLweParams* params){
    const int k = params->k;

    for (int i = 0; i <= k; ++i)
    TorusPolynomial_fft(result->a+i,source->a+i);
    result->current_variance=source->current_variance;
}



EXPORT void tLweFFTClear(TLweSampleFFT* result, const TLweParams* params){
    int k = params->k;

    for (int i = 0; i <= k; ++i) 
    LagrangeHalfCPolynomialClear(&result->a[i]);
    result->current_variance = 0.;
}


EXPORT void tLweFFTAddMulRTo(TLweSampleFFT* result, const LagrangeHalfCPolynomial* p, const TLweSampleFFT* sample, const TLweParams* params) {
    const int k = params->k;

    for (int i=0; i<=k; i++)
    LagrangeHalfCPolynomialAddMul(result->a+i,p,sample->a+i);
    //result->current_variance += sample->current_variance; 
    //TODO: how to compute the variance correctly?
}
























/* Extract functions */
EXPORT void tLweExtractLweSampleIndex(LweSample* result, const TLweSample* x, const int index, const LweParams* params,  const TLweParams* rparams) {
    const int N = rparams->N;
    const int k = rparams->k;
    assert(params->n == k*N);

    for (int i=0; i<k; i++) {
      for (int j=0; j<=index; j++)
        result->a[i*N+j] = x->a[i].coefsT[index-j];
      for (int j=index+1; j<N; j++)
        result->a[i*N+j] = -x->a[i].coefsT[N+index-j];
    }
    result->b = x->b->coefsT[index];
}

EXPORT void tLweExtractLweSample(LweSample* result, const TLweSample* x, const LweParams* params,  const TLweParams* rparams) {
    tLweExtractLweSampleIndex(result, x, 0, params, rparams);
}

//extractions Ring Lwe -> Lwe
EXPORT void tLweExtractKey(LweKey* result, const TLweKey* key) //sans doute un param supplémentaire
{
    const int N = key->params->N;
    const int k = key->params->k;
    assert(result->params->n == k*N);
    for (int i=0; i<k; i++) {
    for (int j=0; j<N; j++)
        result->key[i*N+j]=key->key[i].coefs[j];
    }
}
