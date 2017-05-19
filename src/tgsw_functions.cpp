#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>

using namespace std;



typedef int32_t Torus32; //avant uint32_t
//typedef int64_t Torus64; //avant uint64_t


TGswParams::TGswParams(int l, int Bgbit, const TLweParams* tlwe_params):
    l(l),
    Bgbit(Bgbit),
    Bg(1<<Bgbit),
    halfBg(Bg/2),
    maskMod(Bg-1),
    tlwe_params(tlwe_params),
    kpl(int((tlwe_params->k+1)*l))
    {
        h = new Torus32[l];
        for (int i = 0; i < l; ++i) {
            int kk = (32-(i+1)*Bgbit);
            h[i] = 1 << kk; // 1/(Bg^(i+1)) as a Torus32
        }   

        // offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit)) 
        uint32_t temp1 = 0;
        for (int i = 0; i < l; ++i) {
            uint32_t temp0 = 1 << (32-(i+1)*Bgbit);
            temp1 += temp0;
        }
        offset = temp1*halfBg;

    }
TGswParams::~TGswParams() {delete[] h;}

// same key as in TLwe
TGswKey::TGswKey(const TGswParams* params): params(params),tlwe_params(params->tlwe_params),tlwe_key(tlwe_params){
    key = tlwe_key.key;
}
TGswKey::~TGswKey() {}

TGswSampleFFT::TGswSampleFFT(const TGswParams* params, TLweSampleFFT* all_samples_raw): k(params->tlwe_params->k), l(params->l) {
    all_samples = all_samples_raw;
    sample = new TLweSampleFFT*[(k+1)*l];

    for (int p = 0; p < (k+1); ++p)
    sample[p] = all_samples + p*l;
}
TGswSampleFFT::~TGswSampleFFT() {delete[] sample;}











EXPORT void init_TGswSample(TGswSample* obj, const TGswParams* params) {
    const int k = params->tlwe_params->k;
    const int l = params->l;
    TLweSample* all_sample = new_TLweSample_array((k+1)*l,params->tlwe_params); // tous les samples comme un vecteur ligne
    TLweSample** bloc_sample = new TLweSample*[k+1]; // blocs horizontaux (l lignes) de la matrice TGsw

    for (int p = 0; p < k+1; ++p)
        bloc_sample[p] = all_sample + p*l;

    new(obj) TGswSample(all_sample, bloc_sample, k,l);
}

EXPORT void destroy_TGswSample(TGswSample* obj) {
    const int k = obj->k;
    const int l = obj->l;
    delete_TLweSample_array((k+1)*l,obj->all_sample);
    delete[] obj->bloc_sample;
    obj->~TGswSample();
}


EXPORT void tGswKeyGen(TGswKey* result){
    tLweKeyGen(&result->tlwe_key);
}



EXPORT void tGswClear(TGswSample* result, const TGswParams* params){
    const int kpl = params->kpl;

    for (int p = 0; p < kpl; ++p)
        tLweClear(&result->all_sample[p], params->tlwe_params);
}



EXPORT void tGswAddMuH(TGswSample* result, const IntPolynomial* message, const TGswParams* params) {
    const int k = params->tlwe_params->k;
    const int N = params->tlwe_params->N;
    const int l = params->l;
    const Torus32* h = params->h;
    const int* mu = message->coefs;

    // compute result += H
    for (int bloc = 0; bloc <= k; ++bloc)
    for (int i=0; i<l; i++) {
        Torus32* target = 
        result->bloc_sample[bloc][i].a[bloc].coefsT;
        const Torus32 hi=h[i];
        for (int j=0; j<N; j++) {
        target[j]+=mu[j]*hi;
        }
    }
}


EXPORT void tGswAddMuIntH(TGswSample* result, const int message, const TGswParams* params)
{
    const int k = params->tlwe_params->k;
    const int l = params->l;
    const Torus32* h = params->h;

    // compute result += H
    for (int bloc = 0; bloc <= k; ++bloc)
    for (int i=0; i<l; i++) 
        result->bloc_sample[bloc][i].a[bloc].coefsT[0]+=message*h[i];
}


EXPORT void tGswEncryptZero(TGswSample* result, double alpha, const TGswKey* key){
    const TLweKey* rlkey = &key->tlwe_key;
    const int kpl = key->params->kpl;
        
    for (int p = 0; p < kpl; ++p) {
        tLweSymEncryptZero(&result->all_sample[p], alpha, rlkey);
    }
}







EXPORT void tGswExternMulToTLwe(TLweSample* accum, const TGswSample* sample, const TGswParams* params){
    const TLweParams* par=params->tlwe_params;
    const int N=par->N;
    const int kpl=params->kpl;
    //TODO: improve this new/delete
    IntPolynomial* dec =new_IntPolynomial_array(kpl,N);

    tGswTLweDecompH(dec,accum,params);
    tLweClear(accum,par);
    for (int i=0; i<kpl;i++) {
        tLweAddMulRTo(accum,&dec[i],&sample->all_sample[i],par);
    }

    delete_IntPolynomial_array(kpl, dec);
}


EXPORT void tGswSymEncrypt(TGswSample* result, const IntPolynomial* message, double alpha, const TGswKey* key){
    tGswEncryptZero(result, alpha, key);
    tGswAddMuH(result, message, key->params);
}



EXPORT void tGswSymEncryptInt(TGswSample* result, const int message, double alpha, const TGswKey* key){
    tGswEncryptZero(result, alpha, key);
    tGswAddMuIntH(result, message, key->params);
}




EXPORT void tGswSymDecrypt(IntPolynomial* result, const TGswSample* sample, const TGswKey* key, const int Msize){
    const TGswParams* params = key->params;
    const TLweParams* rlwe_params = params->tlwe_params;
    const int N = rlwe_params->N;
    const int l = params->l;
    const int k = rlwe_params->k;
    TorusPolynomial* testvec = new_TorusPolynomial(N);
    TorusPolynomial* tmp = new_TorusPolynomial(N);
    IntPolynomial* decomp = new_IntPolynomial_array(l,N);
    
    const Torus32 indic = modSwitchToTorus32(1, Msize);
    torusPolynomialClear(testvec);
    testvec->coefsT[0]=indic;
    tGswTorus32PolynomialDecompH(decomp, testvec, params);

    torusPolynomialClear(testvec);
    for (int i=0; i<l; i++) {
    for (int j=1; j<N; j++) assert(decomp[i].coefs[j]==0);
    tLwePhase(tmp, &sample->bloc_sample[k][i], &key->tlwe_key);
    torusPolynomialAddMulR(testvec, decomp+i, tmp);
    }
    for (int i=0; i<N; i++)
    result->coefs[i]=modSwitchFromTorus32(testvec->coefsT[i], Msize);
    
    delete_TorusPolynomial(testvec);
    delete_TorusPolynomial(tmp);
    delete_IntPolynomial_array(l,decomp);
}

//fonction de decomposition
EXPORT void tGswTLweDecompH(IntPolynomial* result, const TLweSample* sample, const TGswParams* params){
    const int k = params->tlwe_params->k;
    const int l = params->l;

    for (int i = 0; i <= k; ++i) // b=a[k]
        tGswTorus32PolynomialDecompH(result+(i*l), &sample->a[i], params);
}





EXPORT void tGswTorus32PolynomialDecompH(IntPolynomial* result, const TorusPolynomial* sample, const TGswParams* params){
    const int N = params->tlwe_params->N;
    const int l = params->l;
    const int Bgbit = params->Bgbit;
    uint32_t* buf = (uint32_t*) sample->coefsT;
//#define __AVX2__ //(to test)
#ifndef __AVX2__
    const uint32_t maskMod = params->maskMod;
    const int32_t halfBg = params->halfBg;
    const uint32_t offset = params->offset;
#else
    const uint32_t* maskMod_addr = &params->maskMod;
    const int32_t* halfBg_addr = &params->halfBg;
    const uint32_t* offset_addr = &params->offset;
    //const uint32_t offset = params->offset;
    //const uint32_t maskMod = params->maskMod;
    //const int32_t halfBg = params->halfBg;
#endif  

    //First, add offset to everyone
#ifndef __AVX2__
    for (int j = 0; j < N; ++j) buf[j]+=offset;
#else
    {
    const uint32_t* sit = buf;
    const uint32_t* send = buf+N;
    __asm__ __volatile__ (
        "vpbroadcastd (%2),%%ymm0\n"
        "0:\n"
        "vmovdqu (%0),%%ymm3\n"
        "vpaddd %%ymm0,%%ymm3,%%ymm3\n" // add offset
        "vmovdqu %%ymm3,(%0)\n"
        "addq $32,%0\n"
        "cmpq %1,%0\n"
        "jb 0b\n"
        : "=r"(sit),"=r"(send),"=r"(offset_addr)
        :  "0"(sit), "1"(send), "2"(offset_addr)
        : "%ymm0","%ymm3","memory"
        );
    }
#endif  
    
    //then, do the decomposition (in parallel)
    for (int p = 0; p < l; ++p)
    {
    const int decal = (32-(p+1)*Bgbit);
#ifndef __AVX2__
    int32_t* res_p = result[p].coefs;
    for (int j = 0; j < N; ++j)
    {
        uint32_t temp1 = (buf[j] >> decal) & maskMod; 
        res_p[j] = temp1 - halfBg;
    }
#else
    int32_t* dst = result[p].coefs;
    const uint32_t* sit = buf;
    const uint32_t* send = buf+N;
    const int32_t* decal_addr = &decal;
    __asm__ __volatile__ (
        "vpbroadcastd (%4),%%ymm0\n"
        "vpbroadcastd (%5),%%ymm1\n"
        "vmovd (%3),%%xmm2\n"
        "1:\n"
        "vmovdqu (%1),%%ymm3\n"
        "VPSRLD %%xmm2,%%ymm3,%%ymm3\n" // shift by decal
        "VPAND %%ymm1,%%ymm3,%%ymm3\n"  // and maskMod
        "VPSUBD %%ymm0,%%ymm3,%%ymm3\n" // sub halfBg
        "vmovdqu %%ymm3,(%0)\n"
        "addq $32,%0\n"
        "addq $32,%1\n"
        "cmpq %2,%1\n"
        "jb 1b\n"
        : "=r"(dst),"=r"(sit),"=r"(send),"=r"(decal_addr),"=r"(halfBg_addr),"=r"(maskMod_addr)
        :  "0"(dst), "1"(sit), "2"(send), "3"(decal_addr), "4"(halfBg_addr) ,"5"(maskMod_addr)
        : "%ymm0","%ymm1","%ymm2","%ymm3","memory"
        );
    /* // verify that the assembly block was ok
    int32_t* res_p = result[p].coefs;
    for (int j = 0; j < N; ++j)
    {
        uint32_t temp1 = (buf[j] >> decal) & maskMod; 
        if (res_p[j] != int32_t(temp1 - halfBg)) {
        fprintf(stderr, "j=%d,buf[j]=%u,decal=%u,mask=%u,halfbg=%d,res_p[j]=%d\n",j,buf[j],decal,maskMod,halfBg,res_p[j]);
        abort();
        }
    }*/

#endif  
    }

    //finally, remove offset to everyone
#ifndef __AVX2__
    for (int j = 0; j < N; ++j) buf[j]-=offset;
#else
    {
    const uint32_t* sit = buf;
    const uint32_t* send = buf+N;
    __asm__ __volatile__ (
        "vpbroadcastd (%2),%%ymm0\n"
        "2:\n"
        "vmovdqu (%0),%%ymm3\n"
        "vpsubd %%ymm0,%%ymm3,%%ymm3\n" // add offset
        "vmovdqu %%ymm3,(%0)\n"
        "addq $32,%0\n"
        "cmpq %1,%0\n"
        "jb 2b\n"
        "vzeroall\n"
        : "=r"(sit),"=r"(send),"=r"(offset_addr)
        :  "0"(sit), "1"(send), "2"(offset_addr)
        : "%ymm0","%ymm3","memory"
        );
    }
#endif  
}




EXPORT void tGswExternProduct(TLweSample* result, const TGswSample* a, const TLweSample* b, const TGswParams* params){
    const TLweParams* parlwe = params->tlwe_params;
    const int N = parlwe->N;
    const int kpl = params->kpl;
    IntPolynomial* dec = new_IntPolynomial_array(kpl,N);

    tGswTLweDecompH(dec, b, params);

    tLweClear(result, parlwe);
    for (int i = 0; i < kpl; i++) 
    tLweAddMulRTo(result, &dec[i], &a->all_sample[i], parlwe);

    result->current_variance += b->current_variance; //todo + the error term?

    delete_IntPolynomial_array(kpl, dec);
}














/* TGSW fft functions --> ILA: de quoi on a vraiment besoin? */

EXPORT void init_TGswSampleFFT(TGswSampleFFT* obj, const TGswParams* params) {
    const int k = params->tlwe_params->k;
    const int l = params->l;
    TLweSampleFFT* all_samples = new_TLweSampleFFT_array((k+1)*l,params->tlwe_params);
    new(obj) TGswSampleFFT(params, all_samples);
}

EXPORT void destroy_TGswSampleFFT(TGswSampleFFT* obj) {
    int k = obj->k;
    int l = obj->l;
    delete_TLweSampleFFT_array((k+1)*l,obj->all_samples);
    obj->~TGswSampleFFT();
}


EXPORT void tGswToFFTConvert(TGswSampleFFT* result, const TGswSample* source, const TGswParams* params) {
    const int kpl = params->kpl;
    
    for (int p=0; p<kpl; p++)
    tLweToFFTConvert(result->all_samples+p, source->all_sample+p, params->tlwe_params);
}

EXPORT void tGswFromFFTConvert(TGswSample* result, const TGswSampleFFT* source, const TGswParams* params){
    const int kpl = params->kpl;
    
    for (int p=0; p<kpl; p++)
    tLweFromFFTConvert(result->all_sample+p, source->all_samples+p, params->tlwe_params);
}



EXPORT void tGswFFTAddH(TGswSampleFFT* result, const TGswParams* params) {
    const int k = params->tlwe_params->k;
    const int l = params->l;

    for (int j=0; j<l; j++) {
        Torus32 hj = params->h[j];
        for (int i=0; i<=k; i++)
       LagrangeHalfCPolynomialAddTorusConstant(&result->sample[i][j].a[i],hj); 
    }

}

EXPORT void tGswFFTClear(TGswSampleFFT* result, const TGswParams* params) {
    const int kpl = params->kpl;

    for (int p=0; p<kpl; p++)
    tLweFFTClear(result->all_samples+p, params->tlwe_params);
}    

EXPORT void tGswFFTExternMulToTLwe(TLweSample* accum, const TGswSampleFFT* gsw, const TGswParams* params) {
    const TLweParams* tlwe_params=params->tlwe_params;
    const int k = tlwe_params->k;
    const int l = params->l;
    const int kpl = params->kpl;
    const int N = tlwe_params->N;
    //TODO attention, improve these new/delete...
    IntPolynomial* deca = new_IntPolynomial_array(kpl,N); //decomposed accumulator 
    LagrangeHalfCPolynomial* decaFFT=new_LagrangeHalfCPolynomial_array(kpl,N); //fft version
    TLweSampleFFT* tmpa = new_TLweSampleFFT(tlwe_params);

    for (int i=0; i<=k; i++)
    tGswTorus32PolynomialDecompH(deca+i*l,accum->a+i, params);
    for (int p=0; p<kpl; p++)
    IntPolynomial_ifft(decaFFT+p,deca+p);

    tLweFFTClear(tmpa, tlwe_params);
    for (int p=0; p<kpl; p++) {
    tLweFFTAddMulRTo(tmpa, decaFFT+p, gsw->all_samples+p, tlwe_params);
    }
    tLweFromFFTConvert(accum, tmpa, tlwe_params);

    delete_TLweSampleFFT(tmpa);
    delete_LagrangeHalfCPolynomial_array(kpl,decaFFT);
    delete_IntPolynomial_array(kpl,deca);
}



