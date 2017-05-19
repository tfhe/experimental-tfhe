#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>

using namespace std;


typedef int32_t Torus32; //avant uint32_t
//typedef int64_t Torus64; //avant uint64_t




/* LWE */
struct LweParams(int n, double alpha_min, double alpha_max) {} n(n), alpha_min(alpha_min), alpha_max(alpha_max);
LweParams::~LweParams() {}

LweSample::LweSample(const LweParams* params) {
    this->a = new Torus32[params->n];
    this->b = 0;
    this->current_variance = 0.;
}
LweSample::~LweSample() {delete[] a;}

LweKey::LweKey(const LweParams* params) {
    this->params = params;
    this->key = new int[params->n];
}
LweKey::~LweKey() {delete[] key;}



void lweKeyGen(LweKey* result) {
  const int n = result->params->n;
  uniform_int_distribution<int> distribution(0,1);

  for (int i=0; i<n; i++) 
    result->key[i]=distribution(generator);
}

void lweSymEncrypt(LweSample* result, Torus32 message, double alpha, const LweKey* key){
    const int n = key->params->n;

    result->b = gaussian32(message, alpha); 
    for (int i = 0; i < n; ++i)
    {
        result->a[i] = uniformTorus32_distrib(generator);
        result->b += result->a[i]*key->key[i];
    }

    result->current_variance = alpha*alpha;
}

Torus32 lwePhase(const LweSample* sample, const LweKey* key){
    const int n = key->params->n;
    Torus32 axs = 0;
    const Torus32 *__restrict a = sample->a;
    const int * __restrict k = key->key;

    for (int i = 0; i < n; ++i) 
       axs += a[i]*k[i]; 
    return sample->b - axs;
}

Torus32 lweSymDecrypt(const LweSample* sample, const LweKey* key, const int Msize){
    Torus32 phi;

    phi = lwePhase(sample, key);
    return approxPhase(phi, Msize);
}


void lweNoiselessTrivial(LweSample* result, Torus32 mu, const LweParams* params){
    const int n = params->n;

    for (int i = 0; i < n; ++i) result->a[i] = 0;
    result->b = mu;
    result->current_variance = 0.;
}













/* Key Switching stuff */
LweKeySwitchKey::LweKeySwitchKey(int n, int t, int basebit, const LweParams* out_params, LweSample* ks0_raw){
    this->basebit=basebit;
    this->out_params=out_params; 
    this->n=n;
    this->t=t;
    this->base=1<<basebit;
    this->ks0_raw = ks0_raw;
    ks1_raw = new LweSample*[n*t];
    ks = new LweSample**[n];

    for (int p = 0; p < n*t; ++p)
        ks1_raw[p] = ks0_raw + base*p;
    for (int p = 0; p < n; ++p)
        ks[p] = ks1_raw + t*p;
}
LweKeySwitchKey::~LweKeySwitchKey() {
    delete[] ks1_raw;
    delete[] ks;
}





void lweCreateKeySwitchKey_fromArray(LweSample*** result, 
    const LweKey* out_key, const double out_alpha, 
    const int* in_key, const int n, const int t, const int basebit){
    const int base=1<<basebit;       // base=2 in [CGGI16]

    for(int i=0;i<n;i++) {
        for(int j=0;j<t;j++){
            for(int k=0;k<base;k++){
                Torus32 x=(in_key[i]*k)*(1<<(32-(j+1)*basebit));
                lweSymEncrypt(&result[i][j][k],x,out_alpha,out_key);
            }
        }
    }
}


void lweKeySwitchTranslate_fromArray(LweSample* result, 
    const LweSample*** ks, const LweParams* params, 
    const Torus32* ai, 
    const int n, const int t, const int basebit){
    const int base=1<<basebit;       // base=2 in [CGGI16]
    const int32_t prec_offset=1<<(32-(1+basebit*t)); //precision
    const int mask=base-1;

    for (int i=0;i<n;i++){
        const uint32_t aibar=ai[i]+prec_offset;
        for (int j=0;j<t;j++){
            const uint32_t aij=(aibar>>(32-(j+1)*basebit)) & mask;
            if(aij != 0) {lweSubTo(result,&ks[i][j][aij],params);}
        }
    }
}



EXPORT void lweCreateKeySwitchKey(LweKeySwitchKey* result, const LweKey* in_key, const LweKey* out_key){
    const int n=result->n;
    const int basebit=result->basebit;
    const int t=result->t;

    lweCreateKeySwitchKey_fromArray(result->ks, out_key, out_key->params->alpha_min, in_key->key, n, t, basebit);
}

EXPORT void lweKeySwitch(LweSample* result, const LweKeySwitchKey* ks, const LweSample* sample){
    const LweParams* params=ks->out_params;
    const int n=ks->n;
    const int basebit=ks->basebit;
    const int t=ks->t;

    lweNoiselessTrivial(result,sample->b,params);
    lweKeySwitchTranslate_fromArray(result, (const LweSample***) ks->ks, params, sample->a, n, t, basebit);
}





EXPORT void init_LweKeySwitchKey(LweKeySwitchKey* obj, int n, int t, int basebit, const LweParams* out_params) {
    const int base=1<<basebit;
    LweSample* ks0_raw = new_LweSample_array(n*t*base, out_params);

    new(obj) LweKeySwitchKey(n,t,basebit,out_params, ks0_raw);
}

EXPORT void destroy_LweKeySwitchKey(LweKeySwitchKey* obj) {
    const int n = obj->n;
    const int t = obj->t;
    const int base = obj->base;
    delete_LweSample_array(n*t*base,obj->ks0_raw);

    obj->~LweKeySwitchKey();
}

EXPORT LweKeySwitchKey* alloc_LweKeySwitchKey() {
    return (LweKeySwitchKey*) malloc(sizeof(LweKeySwitchKey));
}
EXPORT LweKeySwitchKey* alloc_LweKeySwitchKey_array(int nbelts) {
    return (LweKeySwitchKey*) malloc(nbelts*sizeof(LweKeySwitchKey));
}

EXPORT void free_LweKeySwitchKey(LweKeySwitchKey* ptr) {
    free(ptr);
}
EXPORT void free_LweKeySwitchKey_array(int nbelts, LweKeySwitchKey* ptr) {
    free(ptr);
}

//initialize the key structure
//(equivalent of the C++ constructor)
EXPORT void init_LweKeySwitchKey_array(int nbelts, LweKeySwitchKey* obj, int n, int t, int basebit, const LweParams* out_params) {
    for (int i=0; i<nbelts; i++) {
    init_LweKeySwitchKey(obj+i, n,t,basebit,out_params);
    }
}

//destroys the LweKeySwitchKey structure
//(equivalent of the C++ destructor)
EXPORT void destroy_LweKeySwitchKey_array(int nbelts, LweKeySwitchKey* obj) {
    for (int i=0; i<nbelts; i++) {
    destroy_LweKeySwitchKey(obj+i);
    }
}
 
EXPORT LweKeySwitchKey* new_LweKeySwitchKey(int n, int t, int basebit, const LweParams* out_params) {
    LweKeySwitchKey* obj = alloc_LweKeySwitchKey();
    init_LweKeySwitchKey(obj, n,t,basebit,out_params);
    return obj;
}
EXPORT LweKeySwitchKey* new_LweKeySwitchKey_array(int nbelts, int n, int t, int basebit, const LweParams* out_params) {
    LweKeySwitchKey* obj = alloc_LweKeySwitchKey_array(nbelts);
    init_LweKeySwitchKey_array(nbelts, obj, n,t,basebit,out_params);
    return obj;
}

EXPORT void delete_LweKeySwitchKey(LweKeySwitchKey* obj) {
    destroy_LweKeySwitchKey(obj);
    free_LweKeySwitchKey(obj);
}
EXPORT void delete_LweKeySwitchKey_array(int nbelts, LweKeySwitchKey* obj) {
    destroy_LweKeySwitchKey_array(nbelts,obj);
    free_LweKeySwitchKey_array(nbelts,obj);
}















/* Bootstrapping stuff */

LweBootstrappingKey::LweBootstrappingKey(const LweParams* in_out_params, 
    const TGswParams* bk_params,
    const TLweParams* accum_params,
    const LweParams* extract_params,
    TGswSample* bk,
    LweKeySwitchKey* ks): in_out_params(in_out_params), 
    bk_params(bk_params),
    accum_params(accum_params),
    extract_params(extract_params),
    bk(bk), ks(ks) {}
LweBootstrappingKey::~LweBootstrappingKey() {}


LweBootstrappingKeyFFT::LweBootstrappingKeyFFT(const LweParams* in_out_params, 
    const TGswParams* bk_params,
    const TLweParams* accum_params,
    const LweParams* extract_params, 
    const TGswSampleFFT* bkFFT,
    const LweKeySwitchKey* ks): in_out_params(in_out_params), 
    bk_params(bk_params), 
    accum_params(accum_params), 
    extract_params(extract_params),
    bkFFT(bkFFT), ks(ks) {}
LweBootstrappingKeyFFT::~LweBootstrappingKeyFFT() {}
 



EXPORT void init_LweBootstrappingKeyFFT(LweBootstrappingKeyFFT* obj, const LweBootstrappingKey* bk) {
    
    const LweParams* in_out_params = bk->in_out_params;
    const TGswParams* bk_params = bk->bk_params;
    const TLweParams* accum_params = bk_params->tlwe_params;
    const LweParams* extract_params = &accum_params->extracted_lweparams;
    const int n = in_out_params->n;
    const int t = bk->ks->t;
    const int basebit = bk->ks->basebit;
    const int base = bk->ks->base;
    const int N = extract_params->n;

    LweKeySwitchKey* ks = new_LweKeySwitchKey(N, t, basebit, in_out_params);
    // Copy the KeySwitching key
    for(int i=0; i<N; i++) {
        for(int j=0; j<t; j++){
            for(int p=0; p<base; p++) {
                lweCopy(&ks->ks[i][j][p], &bk->ks->ks[i][j][p], in_out_params);
            }
        }
    }

    // Bootstrapping Key FFT 
    TGswSampleFFT* bkFFT = new_TGswSampleFFT_array(n,bk_params);
    for (int i=0; i<n; ++i) {
        tGswToFFTConvert(&bkFFT[i], &bk->bk[i], bk_params);
    }

    new(obj) LweBootstrappingKeyFFT(in_out_params, bk_params, accum_params, extract_params, bkFFT, ks);
}



EXPORT void destroy_LweBootstrappingKeyFFT(LweBootstrappingKeyFFT* obj) {
    delete_LweKeySwitchKey((LweKeySwitchKey*) obj->ks);
    delete_TGswSampleFFT_array(obj->in_out_params->n,(TGswSampleFFT*) obj->bkFFT);
    obj->~LweBootstrappingKeyFFT();
}



void tfhe_MuxRotate_FFT(TLweSample* result, const TLweSample* accum, const TGswSampleFFT* bki, const int barai, const TGswParams* bk_params) {
    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    tLweMulByXaiMinusOne(result, barai, accum, bk_params->tlwe_params);
    tGswFFTExternMulToTLwe(result, bki, bk_params);
    tLweAddTo(result, accum, bk_params->tlwe_params);
}



EXPORT void tfhe_blindRotate_FFT(TLweSample* accum, 
    const TGswSampleFFT* bkFFT, 
    const int* bara, 
    const int n, 
    const TGswParams* bk_params) {

    //TGswSampleFFT* temp = new_TGswSampleFFT(bk_params);
    TLweSample* temp = new_TLweSample(bk_params->tlwe_params);
    TLweSample* temp2 = temp;
    TLweSample* temp3 = accum; 

    for (int i=0; i<n; i++) {
        const int barai=bara[i];
        if (barai==0) continue; //indeed, this is an easy case!
        
        tfhe_MuxRotate_FFT(temp2, temp3, bkFFT+i, barai, bk_params);
        swap(temp2,temp3);
    }
    if (temp3 != accum) {
        tLweCopy(accum, temp3, bk_params->tlwe_params);
    }
    
    delete_TLweSample(temp);
    //delete_TGswSampleFFT(temp);
}




EXPORT void tfhe_blindRotateAndExtract_FFT(LweSample* result, 
    const TorusPolynomial* v, 
    const TGswSampleFFT* bk, 
    const int barb,
    const int* bara,
    const int n,
    const TGswParams* bk_params) {

    const TLweParams* accum_params = bk_params->tlwe_params;
    const LweParams* extract_params = &accum_params->extracted_lweparams;
    const int N = accum_params->N;
    const int _2N = 2*N;

    // Test polynomial 
    TorusPolynomial* testvectbis = new_TorusPolynomial(N);
    // Accumulator
    TLweSample* acc = new_TLweSample(accum_params);

    // testvector = X^{2N-barb}*v
    if (barb!=0) torusPolynomialMulByXai(testvectbis, _2N-barb, v);
    else torusPolynomialCopy(testvectbis, v);
    tLweNoiselessTrivial(acc, testvectbis, accum_params);
    // Blind rotation
    tfhe_blindRotate_FFT(acc, bk, bara, n, bk_params);
    // Extraction
    tLweExtractLweSample(result, acc, extract_params, accum_params);

    delete_TLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}



EXPORT void tfhe_bootstrap_woKS_FFT(LweSample* result, 
    const LweBootstrappingKeyFFT* bk, 
    Torus32 mu, 
    const LweSample* x){

    const TGswParams* bk_params = bk->bk_params;
    const TLweParams* accum_params = bk->accum_params;
    const LweParams* in_params = bk->in_out_params;
    const int N=accum_params->N;
    const int Nx2= 2*N;
    const int n = in_params->n;

    TorusPolynomial* testvect = new_TorusPolynomial(N);
    int* bara = new int[N];
    

    // Modulus switching
    int barb = modSwitchFromTorus32(x->b,Nx2);
    for (int i=0; i<n; i++) {
        bara[i]=modSwitchFromTorus32(x->a[i],Nx2);
    }

    // the initial testvec = [mu,mu,mu,...,mu]
    for (int i=0;i<N;i++) testvect->coefsT[i]=mu;

    // Bootstrapping rotation and extraction
    tfhe_blindRotateAndExtract_FFT(result, testvect, bk->bkFFT, barb, bara, n, bk_params);

    
    delete[] bara;
    delete_TorusPolynomial(testvect);
}


   
EXPORT void tfhe_bootstrap_FFT(LweSample* result, 
    const LweBootstrappingKeyFFT* bk, 
    Torus32 mu, 
    const LweSample* x){

    LweSample* u = new_LweSample(&bk->accum_params->extracted_lweparams);

    tfhe_bootstrap_woKS_FFT(u, bk, mu, x);
    // Key switching
    lweKeySwitch(result, bk->ks, u);

    delete_LweSample(u);
}










EXPORT void init_LweBootstrappingKey(LweBootstrappingKey* obj, int ks_t, int ks_basebit, const LweParams* in_out_params, const TGswParams* bk_params) {
    const TLweParams* accum_params = bk_params->tlwe_params;
    const LweParams* extract_params = &accum_params->extracted_lweparams;
    const int n = in_out_params->n;
    const int N = extract_params->n;
    
    TGswSample* bk=new_TGswSample_array(n,bk_params);
    LweKeySwitchKey* ks=new_LweKeySwitchKey(N, ks_t, ks_basebit, in_out_params);

    new(obj) LweBootstrappingKey(in_out_params, bk_params, accum_params, extract_params, bk, ks);
}
EXPORT void destroy_LweBootstrappingKey(LweBootstrappingKey* obj) {
    delete_LweKeySwitchKey(obj->ks);
    delete_TGswSample_array(obj->in_out_params->n,obj->bk);
    obj->~LweBootstrappingKey();
}



EXPORT void tfhe_createLweBootstrappingKey(
    LweBootstrappingKey* bk, 
    const LweKey* key_in, 
    const TGswKey* rgsw_key) {
    assert(bk->bk_params==rgsw_key->params);
    assert(bk->in_out_params==key_in->params);

    const LweParams* in_out_params = bk->in_out_params; 
    const TGswParams* bk_params = bk->bk_params;
    const TLweParams* accum_params = bk_params->tlwe_params;
    const LweParams* extract_params = &accum_params->extracted_lweparams;

    //LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)
    const TLweKey* accum_key = &rgsw_key->tlwe_key;
    LweKey* extracted_key = new_LweKey(extract_params);
    tLweExtractKey(extracted_key, accum_key);
    lweCreateKeySwitchKey(bk->ks, extracted_key, key_in);
    delete_LweKey(extracted_key);

    //TGswSample* bk; ///< the bootstrapping key (s->s")
    int* kin = key_in->key;
    const double alpha = accum_params->alpha_min;
    const int n = in_out_params->n;
    //const int kpl = bk_params->kpl;
    //const int k = accum_params->k;
    //const int N = accum_params->N;
    //cout << "create the bootstrapping key bk ("  << "  " << n*kpl*(k+1)*N*4 << " bytes)" << endl;
    //cout << "  with noise_stdev: " << alpha << endl;
    for (int i=0; i<n; i++) {
    tGswSymEncryptInt(&bk->bk[i], kin[i], alpha, rgsw_key);
    }
    // renormalize
    // renormalizeBK(bk, key_in, rgsw_key);
}







//allocate memory space for a LweBootstrappingKey
EXPORT LweBootstrappingKey* alloc_LweBootstrappingKey() {
    return (LweBootstrappingKey*) malloc(sizeof(LweBootstrappingKey));
}
EXPORT LweBootstrappingKey* alloc_LweBootstrappingKey_array(int nbelts) {
    return (LweBootstrappingKey*) malloc(nbelts*sizeof(LweBootstrappingKey));
}

//free memory space for a LweKey
EXPORT void free_LweBootstrappingKey(LweBootstrappingKey* ptr) {
    free(ptr);
}
EXPORT void free_LweBootstrappingKey_array(int nbelts, LweBootstrappingKey* ptr) {
    free(ptr);
}

//initialize the key structure
//(equivalent of the C++ constructor)
EXPORT void init_LweBootstrappingKey_array(int nbelts, LweBootstrappingKey* obj, int ks_t, int ks_basebit, const LweParams* in_out_params, const TGswParams* bk_params) {
    for (int i=0; i<nbelts; i++) {
        init_LweBootstrappingKey(obj+i, ks_t, ks_basebit, in_out_params, bk_params);
    }
}

//destroys the LweBootstrappingKey structure
//(equivalent of the C++ destructor)
EXPORT void destroy_LweBootstrappingKey_array(int nbelts, LweBootstrappingKey* obj) {
    for (int i=0; i<nbelts; i++) {
        destroy_LweBootstrappingKey(obj+i);
    }
}
 
//allocates and initialize the LweBootstrappingKey structure
//(equivalent of the C++ new)
EXPORT LweBootstrappingKey* new_LweBootstrappingKey(const int ks_t, const int ks_basebit, const LweParams* in_out_params, const TGswParams* bk_params) {
    LweBootstrappingKey* obj = alloc_LweBootstrappingKey();
    init_LweBootstrappingKey(obj,ks_t,ks_basebit,in_out_params,bk_params);
    return obj;
}
EXPORT LweBootstrappingKey* new_LweBootstrappingKey_array(int nbelts, const int ks_t, const int ks_basebit, const LweParams* in_out_params, const TGswParams* bk_params) {
    LweBootstrappingKey* obj = alloc_LweBootstrappingKey_array(nbelts);
    init_LweBootstrappingKey_array(nbelts,obj,ks_t,ks_basebit,in_out_params,bk_params);
    return obj;
}

//destroys and frees the LweBootstrappingKey structure
//(equivalent of the C++ delete)
EXPORT void delete_LweBootstrappingKey(LweBootstrappingKey* obj) {
    destroy_LweBootstrappingKey(obj);
    free_LweBootstrappingKey(obj);
}
EXPORT void delete_LweBootstrappingKey_array(int nbelts, LweBootstrappingKey* obj) {
    destroy_LweBootstrappingKey_array(nbelts,obj);
    free_LweBootstrappingKey_array(nbelts,obj);
}





//allocate memory space for a LweBootstrappingKeyFFT
EXPORT LweBootstrappingKeyFFT* alloc_LweBootstrappingKeyFFT() {
    return (LweBootstrappingKeyFFT*) malloc(sizeof(LweBootstrappingKeyFFT));
}
EXPORT LweBootstrappingKeyFFT* alloc_LweBootstrappingKeyFFT_array(int nbelts) {
    return (LweBootstrappingKeyFFT*) malloc(nbelts*sizeof(LweBootstrappingKeyFFT));
}

//free memory space for a LweKey
EXPORT void free_LweBootstrappingKeyFFT(LweBootstrappingKeyFFT* ptr) {
    free(ptr);
}
EXPORT void free_LweBootstrappingKeyFFT_array(int nbelts, LweBootstrappingKeyFFT* ptr) {
    free(ptr);
}

//initialize the key structure
EXPORT void init_LweBootstrappingKeyFFT_array(int nbelts, LweBootstrappingKeyFFT* obj, const LweBootstrappingKey* bk) {
    for (int i=0; i<nbelts; i++) {
    init_LweBootstrappingKeyFFT(obj+i,bk);
    }
}

EXPORT void destroy_LweBootstrappingKeyFFT_array(int nbelts, LweBootstrappingKeyFFT* obj) {
    for (int i=0; i<nbelts; i++) {
        destroy_LweBootstrappingKeyFFT(obj+i);
    }
}
 
//allocates and initialize the LweBootstrappingKeyFFT structure
//(equivalent of the C++ new)
EXPORT LweBootstrappingKeyFFT* new_LweBootstrappingKeyFFT(const LweBootstrappingKey* bk) {
    LweBootstrappingKeyFFT* obj = alloc_LweBootstrappingKeyFFT();
    init_LweBootstrappingKeyFFT(obj,bk);
    return obj;
}
EXPORT LweBootstrappingKeyFFT* new_LweBootstrappingKeyFFT_array(int nbelts, const LweBootstrappingKey* bk) {
    LweBootstrappingKeyFFT* obj = alloc_LweBootstrappingKeyFFT_array(nbelts);
    init_LweBootstrappingKeyFFT_array(nbelts,obj,bk);
    return obj;
}

//destroys and frees the LweBootstrappingKeyFFT structure
//(equivalent of the C++ delete)
EXPORT void delete_LweBootstrappingKeyFFT(LweBootstrappingKeyFFT* obj) {
    destroy_LweBootstrappingKeyFFT(obj);
    free_LweBootstrappingKeyFFT(obj);
}
EXPORT void delete_LweBootstrappingKeyFFT_array(int nbelts, LweBootstrappingKeyFFT* obj) {
    destroy_LweBootstrappingKeyFFT_array(nbelts,obj);
    free_LweBootstrappingKeyFFT_array(nbelts,obj);
}

























