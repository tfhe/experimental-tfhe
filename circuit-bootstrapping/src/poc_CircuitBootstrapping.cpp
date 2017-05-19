#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "generic_utils.h"
#define k 1

#include "spqlios/lagrangehalfc_impl.h"

using namespace std;

#include "poc_types.h"

Random* global_random = new Random();

#if 0
// Original paper parameters
const int Globals::n_lvl0=500;
const int Globals::n_lvl1=1024; 
const int Globals::n_lvl2=2048; 
const int Globals::bgbit_lvl1=8;
const int Globals::ell_lvl1=4;
const int Globals::bgbit_lvl2=9;
const int Globals::ell_lvl2=6;
const double Globals::bkstdev_lvl2=pow(2.,-50);
const double Globals::ksstdev_lvl10=pow(2.,-15);
const int Globals::kslength_lvl10=15;
const int Globals::ksbasebit_lvl10=1;
const double Globals::ksstdev_lvl21=pow(2.,-31);
const int Globals::kslength_lvl21=32;
const int Globals::ksbasebit_lvl21=1;
#endif

#if 0 //180 à 210ms
const int Globals::n_lvl0=500;
const int Globals::n_lvl1=1024; 
const int Globals::n_lvl2=2048; 
const int Globals::bgbit_lvl1=8;
const int Globals::ell_lvl1=2;
const int Globals::bgbit_lvl2=9;
const int Globals::ell_lvl2=6;
const double Globals::bkstdev_lvl2=pow(2.,-45);
const double Globals::ksstdev_lvl10=pow(2.,-14);
const int Globals::kslength_lvl10=11;
const int Globals::ksbasebit_lvl10=1;
const double Globals::ksstdev_lvl21=pow(2.,-31);
const int Globals::kslength_lvl21=16; // 32
const int Globals::ksbasebit_lvl21=2; // 1
#endif

#if 0 //155 à 181ms
const int Globals::n_lvl0=500;
const int Globals::n_lvl1=1024; 
const int Globals::n_lvl2=2048; 
const int Globals::bgbit_lvl1=8;
const int Globals::ell_lvl1=2;
const int Globals::bgbit_lvl2=9;
const int Globals::ell_lvl2=4;
const double Globals::bkstdev_lvl2=pow(2.,-45);
const double Globals::ksstdev_lvl10=pow(2.,-14);
const int Globals::kslength_lvl10=6;
const int Globals::ksbasebit_lvl10=2;
const double Globals::ksstdev_lvl21=pow(2.,-31);
const int Globals::kslength_lvl21=16; 
const int Globals::ksbasebit_lvl21=2;
#endif

#if 1 //144 to ???ms
const int Globals::n_lvl0=500;
const int Globals::n_lvl1=1024; 
const int Globals::n_lvl2=2048; 
const int Globals::bgbit_lvl1=8;
const int Globals::ell_lvl1=2;
const int Globals::bgbit_lvl2=9;
const int Globals::ell_lvl2=4;
const double Globals::bkstdev_lvl2=pow(2.,-44);
const double Globals::ksstdev_lvl10=pow(2.,-14);
const int Globals::kslength_lvl10=6;
const int Globals::ksbasebit_lvl10=2;
const double Globals::ksstdev_lvl21=pow(2.,-31);
const int Globals::kslength_lvl21=10; 
const int Globals::ksbasebit_lvl21=3;
#endif


void lwe32Encrypt_lvl0(LweSample32* cipher, const Torus32 mess, const double stdev, const Globals* env){
    const int n = env->n_lvl0;

    *cipher->b = random_gaussian32(mess, stdev);
    for (int i = 0; i < n; ++i) {
        cipher->a[i] = random_int32();
        *cipher->b += cipher->a[i]*env->key_lvl0[i];
    }
}

void lwe32Encrypt_lvl1(LweSample32* cipher, const Torus32 mess, const double stdev, const Globals* env){
    const int n = env->n_lvl1;

    *cipher->b = random_gaussian32(mess, stdev);
    for (int i = 0; i < n; ++i) {
        cipher->a[i] = random_int32();
        *cipher->b += cipher->a[i]*env->key_lvl1[i];
    }
}

Torus32 lwe32Phase_lvl0(const LweSample32* cipher, const Globals* env) {
    const int n = env->n_lvl0;
    Torus32 res = *cipher->b;
    for (int i = 0; i < n; ++i) {
        res -= cipher->a[i]*env->key_lvl0[i];
    }
    return res;
}

Torus32 lwe32Phase_lvl1(const LweSample32* cipher, const Globals* env) {
    const int n = env->n_lvl1;
    Torus32 res = *cipher->b;
    for (int i = 0; i < n; ++i) {
        res -= cipher->a[i]*env->key_lvl1[i];
    }
    return res;
}


Torus64 lwe64Phase_lvl2(const LweSample64* cipher, const Globals* env) {
    const int n = env->n_lvl2;
    Torus64 res = *cipher->b;
    for (int i = 0; i < n; ++i) {
        res -= cipher->a[i]*env->key_lvl2[i];
    }
    return res;
}


void torus32PolynomialMultKaratsuba_lvl1(Torus32Polynomial* result, const IntPolynomial* poly1, const Torus32Polynomial* poly2, const Globals* env);
void torus32PolynomialMultAddKaratsuba_lvl1(Torus32Polynomial* result, const IntPolynomial* poly1, const Torus32Polynomial* poly2, const Globals* env);
void torus64PolynomialMultKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomial* poly1, const Torus64Polynomial* poly2, const Globals* env);
void torus64PolynomialMultAddKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomial* poly1, const Torus64Polynomial* poly2, const Globals* env);


void tLwe32EncryptZero_lvl1(TLweSample32* cipher, const double stdev, const Globals* env){
    const int N = env->N_lvl1;

    for (int j = 0; j < N; ++j) cipher->b->coefs[j] = random_gaussian32(0, stdev);

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < N; ++j) cipher->a[i].coefs[j] = random_int32();
    }
    for (int i = 0; i < k; ++i) torus32PolynomialMultAddKaratsuba_lvl1(cipher->b, &env->Key_lvl1[i], &cipher->a[i], env);
}


void tLwe32Phase_lvl1(Torus32Polynomial* phase, const TLweSample32* cipher, const Globals* env){
    const int N = env->N_lvl1;

    //since we only have AddMult, we compute the opposite of the phase
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -cipher->b->coefs[j];
    }

    for (int i = 0; i < k; ++i) {
        torus32PolynomialMultAddKaratsuba_lvl1(phase, &env->Key_lvl1[i], &cipher->a[i], env);
    }

    //and we negate the result
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -phase->coefs[j];
    }
}

void tLwe64Phase_lvl2(Torus64Polynomial* phase, const TLweSample64* cipher, const Globals* env){
    const int N = env->N_lvl2;

    //since we only have AddMult, we compute the opposite of the phase
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -cipher->b->coefs[j];
    }

    for (int i = 0; i < k; ++i) {
        torus64PolynomialMultAddKaratsuba_lvl2(phase, &env->Key_lvl2[i], &cipher->a[i], env);
    }

    //and we negate the result
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -phase->coefs[j];
    }
}

void tLwe64Encrypt_lvl2(TLweSample64* cipher, const Torus64 mess, const double stdev, const Globals* env){
    const int N = env->N_lvl2;

    cipher->b->coefs[0] = random_gaussian64(mess, stdev);
    for (int j = 1; j < N; ++j) cipher->b->coefs[j] = random_gaussian64(0, stdev);

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < N; ++j) cipher->a[i].coefs[j] = random_int64();
    }
    for (int i = 0; i < k; ++i) torus64PolynomialMultAddKaratsuba_lvl2(cipher->b, &env->Key_lvl2[i], &cipher->a[i], env);

/*
    Torus64Polynomial* phase = new Torus64Polynomial(N);
    tLwe64Phase_lvl2(phase,cipher,env);
    printf("Phase: ");
    for (int i=0; i<N; i++) {
        printf("%ld ",phase->coefs[i]);
    }
    printf("\n");
    delete phase;
*/
}


void tGsw64Encrypt_lvl2(TGswSample64* cipher, const int mess, const double stdev, const Globals* env){
    const int l = env->ell_lvl2;
    const int Bgbit = env->bgbit_lvl2;

    for (int bloc = 0; bloc <= k; ++bloc) {
        for (int i = 0; i < l; ++i) {
            // encryption of 0
            tLwe64Encrypt_lvl2(&cipher->samples[bloc][i], 0, stdev, env);
            // add mess*h[i]
            cipher->samples[bloc][i].a[bloc].coefs[0] += mess * (UINT64_C(1) << (64-(i+1)*Bgbit));
        }
    }    
}






/////////////////////////////////////////////////////////////////
//     FFT SPECIFIC SECTION                                    //
/////////////////////////////////////////////////////////////////

LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int nbelts, int N) {
    return new_array1<LagrangeHalfCPolynomial>(nbelts,N);
}

void delete_LagrangeHalfCPolynomial_array(int nbelts, LagrangeHalfCPolynomial* data) {
    delete_array1<LagrangeHalfCPolynomial>(data);
}


#ifdef USE_FFT
void IntPolynomial_ifft_lvl2(LagrangeHalfCPolynomial* result, const IntPolynomial* source, const Globals* env) {
    assert(env->N_lvl2==2048);
    fftp2048.execute_reverse_int(result->values, source->coefs);
} 


void LagrangeHalfCPolynomialClear_lvl2(LagrangeHalfCPolynomial* result, const Globals* env) {
    const int N = env->N_lvl2;    
    for (int i=0; i<N; i++) 
    result->values[i] = 0;
}

void LagrangeHalfCPolynomialAddMul_lvl2(LagrangeHalfCPolynomial* result, const LagrangeHalfCPolynomial* a, const LagrangeHalfCPolynomial* b, const Globals* env) {
    const int Ns2 = env->N_lvl2/2;
    /*    
    for (int i=0; i<Ns2; i++) {
        double ra = a->values[i];
        double ia = a->values[Ns2+i];
        double rb = b->values[i];
        double ib = b->values[Ns2+i];
        result->values[i] += ra*rb-ia*ib;
        result->values[i+Ns2] += ra*ib+ia*rb;
    }
    */
    LagrangeHalfCPolynomialAddMulASM(result->values, a->values, b->values, Ns2);
}

void TorusPolynomial64_fft_lvl2(Torus64Polynomial* result, const LagrangeHalfCPolynomial* source, const Globals* env) {
    assert(env->N_lvl2==2048);
    fftp2048.execute_direct_torus64(result->coefs, source->values);
} 

void TorusPolynomial64_ifft_lvl2(LagrangeHalfCPolynomial* result, const Torus64Polynomial* source, const Globals* env) {
    assert(env->N_lvl2==2048);
    fftp2048.execute_reverse_torus64(result->values, source->coefs);
}

#else
//these are fake and slow versions of the FFT, that use Karatsuba instead
void IntPolynomial_ifft_lvl2(LagrangeHalfCPolynomial* result, const IntPolynomial* source, const Globals* env) {
    assert(env->N_lvl2==2048);
    result->setIntPoly(source, 2048);
} 


void LagrangeHalfCPolynomialClear_lvl2(LagrangeHalfCPolynomial* result, const Globals* env) {
    assert(env->N_lvl2==2048);
    result->setZeroTorus64Poly(2048);
}

void LagrangeHalfCPolynomialAddMul_lvl2(LagrangeHalfCPolynomial* result, const LagrangeHalfCPolynomial* a, const LagrangeHalfCPolynomial* b, const Globals* env) {
    assert(env->N_lvl2==2048);
    assert(result->torus64Poly!=0);
assert(a->intPoly!=0);
    assert(b->torus64Poly!=0);
    torus64PolynomialMultAddKaratsuba_lvl2(result->torus64Poly, a->intPoly, b->torus64Poly, env); 
}

void TorusPolynomial64_fft_lvl2(Torus64Polynomial* result, const LagrangeHalfCPolynomial* source, const Globals* env) {
    assert(env->N_lvl2==2048);
    assert(source->torus64Poly!=0);
    for (int i=0; i<2048; i++) result->coefs[i]=source->torus64Poly->coefs[i];
} 

void TorusPolynomial64_ifft_lvl2(LagrangeHalfCPolynomial* result, const Torus64Polynomial* source, const Globals* env) {
    assert(env->N_lvl2==2048);
    result->setTorus64Poly(source, 2048);
}
#endif


/////////////////////////////////////////////////////////////////
//     FIN FFT SPECIFIC SECTION                                    //
/////////////////////////////////////////////////////////////////







    
    
   










Globals::Globals(){
        t_lvl0 = kslength_lvl10 * ksbasebit_lvl10; 
        t_lvl1 = kslength_lvl21 * ksbasebit_lvl21;
        N_lvl1 = n_lvl1;
        N_lvl2 = n_lvl2;
        
        // Offset
        torusDecompOffset = 0;
        for (int i = 0; i <= ell_lvl2; ++i) torusDecompOffset |= (UINT64_C(1)<<(63-i*bgbit_lvl2));
        
        // Buffer
        torusDecompBuf = new uint64_t[N_lvl2]; 
        
        //secret keys (and their polynomial interpretation)
        // level 0
        key_lvl0 = new int[n_lvl0];
        for (int i = 0; i < n_lvl0; ++i) key_lvl0[i] = random_bit();
        // level 1
        key_lvl1 = new int[n_lvl1];
        for (int i = 0; i < n_lvl1; ++i) key_lvl1[i] = random_bit();
        Key_lvl1 = new IntPolynomial(N_lvl1);
        for (int i = 0; i < N_lvl1; ++i) Key_lvl1->coefs[i] = key_lvl1[i]; 
        // level 2
        key_lvl2 = new int[n_lvl2+1];
        for (int i = 0; i < n_lvl2; ++i) key_lvl2[i] = random_bit();
        key_lvl2[n_lvl2] = -1; // coef extended key (phase in 1 loop, privKS)
        Key_lvl2 = new IntPolynomial(N_lvl2);
        for (int i = 0; i < N_lvl2; ++i) Key_lvl2->coefs[i] = key_lvl2[i]; 
        
        //cloud keys
        //PreKS: LweSamples at level 0 of the key key_lvl1  encrypted with key_lvl0
        printf("generating preKS...\n");
        int ksbase_lvl10 = 1 << ksbasebit_lvl10;
        preKS = new_array3<LweSample32>(n_lvl1, kslength_lvl10, ksbase_lvl10, n_lvl0); // preKS[n_lvl1][kslength_lvl10][ksbase_lvl10]
        for (int i = 0; i < n_lvl1; ++i) {
            for (int j = 0; j < kslength_lvl10; ++j) {
                for (int u = 0; u < ksbase_lvl10; ++u) {
                    Torus32 messPreKS = (key_lvl1[i] << (32-(j+1)*ksbasebit_lvl10)) * u;
                    lwe32Encrypt_lvl0(&preKS[i][j][u], messPreKS, ksstdev_lvl10, this);
                } 
            }    
        }

#if 1 //TODO
        //circuit bootstrapping keys
        printf("generating bk...\n");
        bk = new_array1<TGswSample64>(n_lvl0, ell_lvl2, N_lvl2); // bk[n_lvl0]
        for (int i = 0; i < n_lvl0; ++i) {
            tGsw64Encrypt_lvl2(&bk[i], key_lvl0[i], bkstdev_lvl2, this);
        }


        printf("generating bkFFT...\n");
        bkFFT = new_array1<TGswSampleFFT>(n_lvl0, ell_lvl2, N_lvl2); // bkFFT[n_lvl0]
        for (int i = 0; i < n_lvl0; ++i) {
            for (int j = 0; j < ((k+1)*ell_lvl2); ++j){
                for (int p = 0; p <= k; ++p) {
                    TorusPolynomial64_ifft_lvl2(&bkFFT[i].allsamples[j].a[p], &bk[i].allsamples[j].a[p], this);
                }
            }
        }


        //privKS
        printf("generating privKS...\n");
        int ksbase_lvl21 = 1 << ksbasebit_lvl21;
        privKS = new_array4<TLweSample32>(k+1, n_lvl2+1, kslength_lvl21, ksbase_lvl21, N_lvl1);
        for (int z = 0; z <= k; ++z) {
            for (int i = 0; i <= n_lvl2; ++i) {
                for (int j = 0; j < kslength_lvl21; ++j) {
                    for (int u = 0; u < ksbase_lvl21; ++u) {
                        Torus32 messPrivKS = (key_lvl2[i] << (32-(j+1)*ksbasebit_lvl21)) * u;
                        tLwe32EncryptZero_lvl1(&privKS[z][i][j][u], ksstdev_lvl21, this);
                        privKS[z][i][j][u].a[z].coefs[0] += messPrivKS;
                    }
                }
            }
        }
#endif


}

// Globals::~Globals(){}








//PreKeySwitch and ModRescale
// result is at level 0 while x is at level 1
// the KS key (preKS) is at level 0
void preKeySwitch(LweSample32* result, const LweSample32* x, const Globals* env) {
    const int n_lvl0 = env->n_lvl0;
    const int n_lvl1 = env->n_lvl1;
    const int kslen = env->kslength_lvl10;
    const int basebit10 = env->ksbasebit_lvl10;
    const int base10 = 1<<basebit10;       // base=2 in [CGGI16]
    const int mask = base10 - 1;
    const int32_t prec_offset = 1<<(32-(1+basebit10*kslen)); //precision
    
    // lwe noiseless trivial (0,b)
    *result->b = *x->b;
    for (int i = 0; i < n_lvl0; ++i) {
        result->a[i] = 0;
    }

    // key switching
    for (int i = 0; i < n_lvl1; ++i) {
        const uint32_t aibar = x->a[i] + prec_offset;

        for (int j = 0; j < kslen; ++j) {
            const uint32_t aij = (aibar>>(32-(j+1)*basebit10)) & mask;
            
            if(aij != 0) {
                for (int h = 0; h <= n_lvl0; ++h) result->a[h] -= env->preKS[i][j][aij].a[h];
            }
        }
    }
    
}




// mod rescale the sample x between 0 and 2.n_lvl2
// x is at level 0
void preModSwitch(int* result, const LweSample32* x, const Globals* env) {
    const int n_lvl0 = env->n_lvl0;
    const int N_lvl2 = env->n_lvl2; // N_lvl2 = n_lvl2
    const int _2N = 2*N_lvl2;
    uint64_t interv = ((UINT64_C(1)<<63)/_2N)*2; // width of each intervall
    uint64_t half_interval = interv/2; // begin of the first intervall

    // Mod Switching (as in modSwitchFromTorus32)
    for (int i = 0; i <= n_lvl0; ++i){
        uint64_t temp = (uint64_t(x->a[i])<<32) + half_interval; // RIVEDI
        result[i] = temp/interv;
    }
}







void tGswTorus64PolynomialDecompH(IntPolynomial* result, const Torus64Polynomial* sample, const Globals* env){
    const int N = env->N_lvl2;
    const int l = env->ell_lvl2;
    const int Bgbit = env->bgbit_lvl2;
    // peut etre tout cela dans le env
    const uint64_t Bg = UINT64_C(1)<<Bgbit; 
    const uint64_t mask = Bg-1;
    const int32_t halfBg = Bg/2; 
    uint64_t* buf = env->torusDecompBuf;
    const uint64_t offset = env->torusDecompOffset;

    //First, add offset to everyone
    for (int j = 0; j < N; ++j) buf[j]=sample->coefs[j]+offset;
    
    //then, do the decomposition (in parallel)
    for (int p = 0; p < l; ++p) {
        const int decal = (64-(p+1)*Bgbit);
        int32_t* res_p = result[p].coefs; // res is a int (ok 32)
        for (int j = 0; j < N; ++j) {
            uint32_t temp1 = (buf[j] >> decal) & mask; 
            res_p[j] = temp1 - halfBg;
        }
    }
}



void tGsw64DecompH(IntPolynomial* result, const TLweSample64* sample, const Globals* env){
    const int l = env->ell_lvl2;
    for (int i = 0; i <= k; ++i) tGswTorus64PolynomialDecompH(result+(i*l), &sample->a[i], env);
}





//circuit bootstrap without keyswitch (outcome -mu,+mu depending on the sign of phibar)
// result is at level 2 while abar is the preModSwitch of sample x at level 0
void circuitBootstrapWoKS(LweSample64* result, const Torus64 mu, const int* abar, const Globals* env) {
    const int N_lvl2 = env->N_lvl2; 
    const int n_lvl0 = env->n_lvl0;
    const int l = env->ell_lvl2;
    const int N2 = N_lvl2/2;
    const int _2l = 2*l;
    const Torus64 mu2 = mu/2;
    Torus64Polynomial* testvecttemp = new Torus64Polynomial(N_lvl2); 
    Torus64Polynomial* testvect = new Torus64Polynomial(N_lvl2); 
#ifdef PARANO
    Torus64Polynomial* phase = new Torus64Polynomial(N_lvl2); //TODO 
#endif
    const int bbar = abar[n_lvl0]; 
    TLweSample64* acc1 = new TLweSample64(N_lvl2);  
    TLweSample64* acc2 = new TLweSample64(N_lvl2); 
    TLweSample64* acc = new TLweSample64(N_lvl2); 
    TLweSampleFFT* accFFT = new TLweSampleFFT(N_lvl2);

    TGswSampleFFT* bkFFT = env->bkFFT; // bkFFT[n_lvl0]
    
    
    // Test Vector = (1+X+...+X^{N-1})*X^{N/2}*mu2
    for (int j = 0; j < N2; ++j) testvecttemp->coefs[j] = -mu2;
    for (int j = N2; j < N_lvl2; ++j) testvecttemp->coefs[j] = mu2;
    // Test Vector *= X^{bbar}
    if (bbar < N_lvl2) {
        for (int j = 0; j < bbar; ++j) testvect->coefs[j] = -testvecttemp->coefs[j-bbar+N_lvl2];
        for (int j = bbar; j < N_lvl2; ++j) testvect->coefs[j] = testvecttemp->coefs[j-bbar];
    } else {
        int bbis = bbar-N_lvl2;
        for (int j = 0; j < bbis; ++j) testvect->coefs[j] = testvecttemp->coefs[j-bbis+N_lvl2];
        for (int j = bbis; j < N_lvl2; ++j) testvect->coefs[j] = -testvecttemp->coefs[j-bbis];
    }

    // acc = encrypt of testvect as a noiseless trivial TLweSample
    for (int j = 0; j < N_lvl2; ++j) {
        acc->a[0].coefs[j] = 0; // (k=1)
        acc->a[1].coefs[j] = testvect->coefs[j]; // b = &a[k]
    }
    
#ifdef PARANO    //TODO: verif
    tLwe64Phase_lvl2(phase, acc, env);
    printf("i: ");
    for (int j=0; j<N_lvl2; j++) printf("%d|%ld ",(phase->coefs[j]>0?1:-1),abs(phase->coefs[j])-mu2);
    printf("\n");
#endif

    IntPolynomial* decomp = new_array1<IntPolynomial>(_2l,N_lvl2); // new_IntPolynomial_array(_2l,N_lvl2);
    LagrangeHalfCPolynomial* decompFFT = new_array1<LagrangeHalfCPolynomial>(_2l,N_lvl2); // new_LagrangeHalfCPolynomial_array(_2l,N_lvl2);   
    // blind rotation 
    for (int i = 0; i < n_lvl0; ++i) {
        int aibar = abar[i];
        if (aibar == 0) continue;

        // acc = BKi*[(X^aibar-1)*acc1]+acc

        // acc1 = acc  
        for (int q = 0; q <= k; ++q) 
            for (int j = 0; j < N_lvl2; ++j) acc1->a[q].coefs[j] = acc->a[q].coefs[j]; 

        // acc2 = (X^aibar-1)*acc1 = acc1*X^aibar - acc1
        for (int q = 0; q <= k; ++q) {
            if (aibar < N_lvl2) {
                for (int j = 0; j < aibar; ++j) acc2->a[q].coefs[j] = acc1->a[q].coefs[j-aibar+N_lvl2] - acc1->a[q].coefs[j];
                for (int j = aibar; j < N_lvl2; ++j) acc2->a[q].coefs[j] = acc1->a[q].coefs[j-aibar] - acc1->a[q].coefs[j];
            } else {
                for (int j = 0; j < aibar-N_lvl2; ++j) acc2->a[q].coefs[j] = acc1->a[q].coefs[j-aibar+N_lvl2] - acc1->a[q].coefs[j];
                for (int j = aibar-N_lvl2; j < N_lvl2; ++j) acc2->a[q].coefs[j] = -acc1->a[q].coefs[j-aibar] - acc1->a[q].coefs[j];
            }
        }

#ifdef PARANO        //TODO: verif
        tLwe64Phase_lvl2(phase, acc2, env);
        printf("ii=%d: ",i);
        for (int j=0; j<N_lvl2; j++) printf("%d|%ld ",(phase->coefs[j]>0?1:-1),abs(phase->coefs[j])-mu2);
        printf("\n");
#endif


        // acc1 = BKi*acc2 (tGswFFTExternMulToTLwe with bkFFT)
        // decompose acc2
        tGsw64DecompH(decomp, acc2, env); // SEE IT
        // decomp in FFT mode 
        for (int p = 0; p < _2l; ++p) IntPolynomial_ifft_lvl2(decompFFT+p,decomp+p, env); 
        // accFFT initialization 
        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFT->a+q, env);
        // external product FFT
        for (int p = 0; p < _2l; ++p) 
            for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a+q, decompFFT+p, &bkFFT->allsamples[p].a[q], env);
        // conversion from FFT
        for (int q = 0; q <= k; ++q) TorusPolynomial64_fft_lvl2(acc1->a+q,accFFT->a+q, env);

#ifdef PARANO        //TODO: verif 2
        tLwe64Phase_lvl2(phase, acc1, env);
        printf("iii=%d: ",i);
        for (int j=0; j<N_lvl2; j++) printf("%ld ",phase->coefs[j]);
        printf("\n");
#endif


        // acc += acc1
        for (int q = 0; q <= k; ++q)
            for (int j = 0; j < N_lvl2; ++j) acc->a[q].coefs[j] += acc1->a[q].coefs[j]; 

#ifdef PARANO        //TODO: verif 2
        //TODO: verif
        tLwe64Phase_lvl2(phase, acc, env);
        printf("i=%d: ",i);
        for (int j=0; j<N_lvl2; j++) printf("%d|%ld ",(phase->coefs[j]>0?1:-1),abs(phase->coefs[j])-mu2);
        printf("\n");
#endif

    }

    
    // Sample Extraction
    result->a[0] = acc->a[0].coefs[0];
    for (int j=1; j<N_lvl2; j++) result->a[j] = -acc->a[0].coefs[N_lvl2-j];
    *result->b = acc->a[1].coefs[0] + mu2;

        
    delete accFFT;
    delete_array1<LagrangeHalfCPolynomial>(decompFFT);
    delete_array1<IntPolynomial>(decomp);
    delete acc;
    delete acc2;
    delete acc1;
    delete testvect;
    delete testvecttemp;
} 





//computes a sample of Ku.x, where Ku = Key_lvl1 if u==0 and Ku = 1 if u==1 
// result is at level 1 while x is at level 2
void circuitPrivKS(TLweSample32* result, const int u, const LweSample64* x, const Globals* env) {
    const int n_lvl2 = env->n_lvl2;
    const int kslen = env->kslength_lvl21;
    const int N_lvl1 = env->n_lvl1; // N_lvl1 = n_lvl1
    const int basebit21 = env->ksbasebit_lvl21;
    const int base21 = 1<<basebit21;       // base=2 in [CGGI16]
    const int mask = base21 - 1;
    const int64_t prec_offset = UINT64_C(1)<<(64-(1+basebit21*kslen)); //precision ILA: revoir

    // clear result
    for (int i = 0; i <= k ; ++i) {
        for (int j = 0; j < N_lvl1; ++j) {
            result->a[i].coefs[j] = 0;
        }
    }

    // Private Key Switching 
    for (int i = 0; i <= n_lvl2; ++i) {
        const uint64_t aibar = x->a[i] + prec_offset;

        for (int j = 0; j < kslen; ++j) {
            const uint64_t aij = (aibar>>(64-(j+1)*basebit21)) & mask;

            if (aij != 0){
                for (int q = 0; q <= k; ++q) {
                    for (int p = 0; p < N_lvl1; ++p) result->a[q].coefs[p] -= env->privKS[u][i][j][aij].a[q].coefs[p];
                }
            }
        }
    }

} 



// LweBootstrappingKeyCircuit contient 
//  * une KS (LWE to LWE)
//  * une BK (LWE to LWE) et BKfft --> modifier l'algo de bootstrap car ils apparaissent les coefs 1/Bg^w!!!
//  * k+1 KS (LWE to TLWE)


// KS (LWE to LWE)
// LweKeySwitchKey(int n, int t, int basebit, const LweParams* out_params, LweSample* ks0_raw);

// k+1 KS (LWE to TLWE)
/*
LweKeySwitchKeyCircuit::LweKeySwitchKeyCircuit(int k, int n, int t, int basebit, const TLweParams* tlwe_params, TLweSample* ks0_raw){
    this->basebit=basebit;
    this->tlwe_params=tlwe_params; 
    this->k=k;
    this->n=n;
    this->t=t;
    this->base=1<<basebit;
    this->ks0_raw = ks0_raw;
    // peut etre avec TLweSample_array
    ks1_raw = new TLweSample*[n*t*(k+1)];
    ks2_raw = new TLweSample**[n*t];
    ks = new TLweSample***[n];

    for (int p = 0; p < n*t*(k+1); ++p)
        ks1_raw[p] = ks0_raw + base*p;
    for (int p = 0; p < n*t; ++p)
        ks2_raw[p] = ks1_raw + (k+1)*p;
    for (int p = 0; p < n; ++p)
        ks[p] = ks2_raw + t*p;
}

LweKeySwitchKeyCircuit::~LweKeySwitchKeyCircuit() {
    delete[] ks2_raw;
    delete[] ks1_raw;
    delete[] ks;
}
*/


/*
void lweCreateKeySwitchKeyCircuit(LweKeySwitchKeyCircuit* result, const LweKey* in_key, const TLweKey* out_key){
    const int k=result->k;
    const int n=result->n;
    const int t=result->t;
    const int basebit=result->basebit;
    const int base=1<<basebit;

    for (int u = 0; u <= k; ++u) { // "position" in the TLweSample
        for (int i = 0; i < n; ++i) { // bit key
            for (int j = 0; j < t; ++j) { // power of 2
                for (int h = 0; h < base; ++h) {
                    Torus32 x = (in_key[i]*h)*(1<<(32-(j+1)*basebit));
                    tLweSymEncryptZero(&result[u][i][j][h], out_key->params->alpha_min, out_key);
                    result[u][i][j][h]->a[u] += x;
                }
            }
        }
    }

}
*/





/*
LweBootstrappingKeyCircuit::LweBootstrappingKeyCircuit(const LweParams* in_out_params, 
    const TGswParams* bk_params,
    const TLweParams* accum_params,
    const LweParams* extract_params,
    TGswSample* bk,
    LweKeySwitchKey* ks);

LweBootstrappingKeyFFT::LweBootstrappingKeyFFT(const LweParams* in_out_params, 
    const TGswParams* bk_params,
    const TLweParams* accum_params,
    const LweParams* extract_params, 
    const TGswSampleFFT* bkFFT,
    const LweKeySwitchKey* ks);
*/

// tfhe_createLweBootstrappingKey(LweBootstrappingKeyCircuit* bk, const LweKey* key_in, const TGswKey* rgsw_key);




/*
void lweToTLweKeySwitching(TLweSample* result, const int u, const LweKeySwitchKeyCircuit* ksCircuit, const LweSample* sample){
    const TLweParams* params = ksCircuit->tlwe_params;
    const int n = ksCircuit->n;
    const int basebit = ksCircuit->basebit;
    const int t = ksCircuit->t;
    const int32_t prec_offset = 1<<(32-(1+basebit*t)); //precision
    const int base = 1<<basebit;
    const int mask = base-1;

    tLweNoiselessTrivialT(result, sample->b, params);

    for (int i = 0; i < n; i++) {
        const uint32_t aibar = sample->a[i] + prec_offset;
        
        for (int j = 0; j < t; ++j) {
            const uint32_t aij = (aibar>>(32-(j+1)*basebit)) & mask;
            if(aij != 0) {tLweSubTo(result, &ks[u][i][j][aij], params);}
        }
    }

}
*/










void tfhe_CircuitBootstrapFFT(TGswSample32* result, const LweSample32* sample, const Globals* env){
    const int ell1 = env->ell_lvl1;
    const int bgbit1 = env->bgbit_lvl1;
    const int n0 = env->n_lvl0;
    const int n2 = env->n_lvl2;
    

    // Pre Key switch from LweSample level 1 to LweSample level 0
    LweSample32* res_preKS = new LweSample32(n0);
    preKeySwitch(res_preKS, sample, env);

    // Bootstrapping (with preModSwitch) from LweSample level 0 to LweSample level 2
    int* res_preMS = new int[n0+1];  // the +1 is for the b part
    preModSwitch(res_preMS, res_preKS, env);
#ifndef NDEBUG
    printf("prems: ");
    for (int i=0; i<n0+1; i++)
        printf("%d ",res_preMS[i]);
    printf("\n");
#endif

    LweSample64* res_boot = new LweSample64(n2);       
    for (int w = 0; w < ell1; ++w) {
        const Torus64 mu1 = UINT64_C(1)<<(64-(w+1)*bgbit1);
        circuitBootstrapWoKS(res_boot, mu1, res_preMS, env);
#ifndef NDEBUG
        cout << mu1 << " " << lwe64Phase_lvl2(res_boot, env) << endl;
#endif
        // Key Switching to fill the TGswSample from LweSample level 2 to TGswSample level 1
        for (int u = 0; u <= k; ++u) {
            //int index = w*(k+1) + u;
            circuitPrivKS(&result->samples[u][w], u, res_boot, env);
        } 
#ifndef NDEBUG 
        for (int u = 0; u <=k; ++u) {
            Torus32Polynomial* phasePoly = new Torus32Polynomial(env->N_lvl1);
            tLwe32Phase_lvl1(phasePoly, &result->samples[u][w], env);
            cout << u << " " << w << " - ";
            for (int j = 0; j < env->N_lvl1; ++j) {
                cout << phasePoly->coefs[j] << " ";
            }
            cout << endl;
            delete phasePoly;
        }
#endif
    }

    delete res_boot;
    delete res_preMS;
    delete res_preKS;
}



/*
void CMux(TLweSample32* out, const TGswSample32* c, const TLweSample32* in0, const TLweSample32* in1, const Globals* env){}
*/
























/* ***********************************************************************************
********************************** MAIN **********************************************
*********************************************************************************** */





int main(int argc, char** argv)
{
/*
    double valsd = 123456*pow(2,50);
    double* valsp = &valsd;
    uint64_t valsi = *((uint64_t*) valsp);
    static const uint64_t valmask0 = 0x000FFFFFFFFFFFFFul;
    static const uint64_t valmask1 = 0x0010000000000000ul;
    static const uint16_t expmask0 = 0x07FFu;
    uint64_t val = (valsi&valmask0)|valmask1; //mantissa on 53 bits
    uint16_t expo = (valsi>>52)&expmask0; //exponent 11 bits
    // 1023 -> 52th pos -> 0th pos
    // 1075 -> 52th pos -> 52th pos
    int16_t trans = expo-1075;
    uint64_t val2 = trans>0?(val<<trans):(val>>-trans);
    int64_t resi=(valsi>>63)?-val2:val2;
    printf("%lf -> %ld\n",valsd,resi);
    exit(0);
*/
    cout << "generating the keys" << endl;
    Globals* env = new Globals();

/*

    int N2 = 2048;
    Torus64Polynomial* pa = new Torus64Polynomial(N2);
    IntPolynomial* pb = new IntPolynomial(N2);
    Torus64Polynomial* pc = new Torus64Polynomial(N2);
    Torus64Polynomial* pd = new Torus64Polynomial(N2);
    LagrangeHalfCPolynomial* pla = new LagrangeHalfCPolynomial(N2);
    LagrangeHalfCPolynomial* plb = new LagrangeHalfCPolynomial(N2);
    LagrangeHalfCPolynomial* plc = new LagrangeHalfCPolynomial(N2);
    for (int i=0; i<N2; i++) 
        pa->coefs[i]=random_int64();
    for (int i=0; i<N2; i++) 
        pb->coefs[i]=random_int32()>>23;
    //sage output:
    printf("sage: pa = [");
    for (int i=0; i<N2; i++) 
        printf("%ld, ",pa->coefs[i]);
    printf("0]\n");
    printf("sage: pb = [");
    for (int i=0; i<N2; i++) 
        printf("%d, ",pb->coefs[i]);
    printf("0]\n");
    
    //FFT
    TorusPolynomial64_ifft_lvl2(pla, pa, env);
    printf("sage: pla = [");
    for (int i=0; i<N2; i++) 
        printf("%lf, ",pla->values[i]);
    printf("0]\n");
    IntPolynomial_ifft_lvl2(plb, pb, env);
    printf("sage: plb = [");
    for (int i=0; i<N2; i++) 
        printf("%lf, ",plb->values[i]);
    printf("0]\n");
    //test reverse
    TorusPolynomial64_fft_lvl2(pc,pla,env);
    for (int i=0; i<N2; i++) 
        printf("--- %ld %ld || %ld\n", pa->coefs[i],pc->coefs[i],pa->coefs[i]-pc->coefs[i] );
    TorusPolynomial64_fft_lvl2(pc,plb,env);
    for (int i=0; i<N2; i++) 
        printf("xxx %d %ld || %ld\n", pb->coefs[i],pc->coefs[i],pb->coefs[i]-pc->coefs[i] );

    LagrangeHalfCPolynomialClear_lvl2(plc, env);
    LagrangeHalfCPolynomialAddMul_lvl2(plc, pla, plb, env);
    TorusPolynomial64_fft_lvl2(pc, plc, env);
    printf("sage: plc = [");
    for (int i=0; i<N2; i++) 
        printf("%lf, ",plc->values[i]);
    printf("0]\n");

    //Karatsuba
    torus64PolynomialMultKaratsuba_lvl2(pd,pb,pa,env);
    for (int i=0; i<N2; i++) 
        printf("%ld %ld || %ld\n", pd->coefs[i],pc->coefs[i],pd->coefs[i]-pc->coefs[i] );
    exit(0);
   
*/

    LweSample32* in_sample = new LweSample32(env->N_lvl1);

    TGswSample32* out_sample = new TGswSample32(env->ell_lvl1, env->N_lvl1);
    
    int32_t message = int32_t(3)<<29;
    lwe32Encrypt_lvl1(in_sample, message, 0.01, env);
    printf("Input message: %d, phase: %d\n", message, lwe32Phase_lvl1(in_sample,env) );

#ifdef NDEBUG
    static const int NBTRIALS=10;
#else
    static const int NBTRIALS=1;
#endif

    cout << "starting circuit bootstrapping " << endl;
    clock_t begin = clock();
    for (int i = 0; i < NBTRIALS; ++i)
    {
        // Circuit bootstrapping
        tfhe_CircuitBootstrapFFT(out_sample, in_sample, env);
    }
    clock_t end = clock();
    cout << "finished circuit bootstrapping " << endl;
    cout << "total time (microsecs for " << NBTRIALS << " TRIALS)... " << (end-begin) << endl;





    return 0;
}































