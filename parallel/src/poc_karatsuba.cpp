#include "poc_types.h" 

#define k 1

/* ***********************************************************************************
******************************* Karatsuba 32 *****************************************
*********************************************************************************** */

void torus32PolynomialMultNaive_plain_aux(Torus32* __restrict result, const int* __restrict poly1, const Torus32* __restrict poly2, const int N) {
    for (int i = 0; i < N; ++i) {
        Torus32 ri = 0;
        for (int j = 0; j <= i; ++j) {
            ri += poly1[j] * poly2[i-j];
        }
        for (int j = i+1; j < N; ++j) {
            ri -= poly1[j] * poly2[N+i-j];
        }
        result[i] = ri;
    }
}


// A and B of size = size
// R of size = 2*size-1
void Karatsuba32_aux(Torus32* R, const int* A, const Torus32* B, const int size, const char* buf){
    const int h = size / 2;
    const int sm1 = size-1;

    //we stop the karatsuba recursion at h=4, because on my machine,
    //it seems to be optimal
    if (h<=4) {
        torus32PolynomialMultNaive_plain_aux(R, A, B, size);
        return;
    }

    //we split the polynomials in 2
    int* Atemp = (int*) buf;
    buf += h*sizeof(int);
    Torus32* Btemp = (Torus32*) buf;
    buf+= h*sizeof(Torus32);
    Torus32* Rtemp = (Torus32*) buf;
    buf+= size*sizeof(Torus32);
    //Note: in the above line, I have put size instead of sm1 so that buf remains aligned on a power of 2

    for (int i = 0; i < h; ++i) Atemp[i] = A[i] + A[h+i];
    for (int i = 0; i < h; ++i) Btemp[i] = B[i] + B[h+i];

    // Karatsuba recursivly
    Karatsuba32_aux(R, A, B, h, buf); // (R[0],R[2*h-2]), (A[0],A[h-1]), (B[0],B[h-1])
    Karatsuba32_aux(R+size, A+h, B+h, h, buf); // (R[2*h],R[4*h-2]), (A[h],A[2*h-1]), (B[h],B[2*h-1])
    Karatsuba32_aux(Rtemp, Atemp, Btemp, h, buf);
    R[sm1]=0; //this one needs to be set manually
    for (int i = 0; i < sm1; ++i) Rtemp[i] -= (R[i] + R[size+i]);
    for (int i = 0; i < sm1; ++i) R[h+i] += Rtemp[i];
}




// ILA: modif const Globals* env --> const int N1
// poly1, poly2 and result are polynomials mod X^N+1
void torus32PolynomialMultKaratsuba_lvl1(Torus32Polynomial* result, const IntPolynomial* poly1, const Torus32Polynomial* poly2, const int N1){
    //const int N1 = env->N_lvl1;
    Torus32* R = new Torus32[2*N1-1];
    char* buf = new char[16*N1]; //that's large enough to store every tmp variables (2*2*N*4)
    
    // Karatsuba 
    Karatsuba32_aux(R, poly1->coefs, poly2->coefs, N1, buf);

    // reduction mod X^N+1
    for (int i = 0; i < N1-1; ++i) result->coefs[i] = R[i] - R[N1+i];
    result->coefs[N1-1] = R[N1-1];
    
    delete[] R;
    delete[] buf;
}




void torus32PolynomialMultAddKaratsuba_lvl1(Torus32Polynomial* result, const IntPolynomial* poly1, const Torus32Polynomial* poly2, const Globals* env){
    const int N1 = env->N_lvl1;
    Torus32* R = new Torus32[2*N1-1];
    char* buf = new char[16*N1]; //that's large enough to store every tmp variables (2*2*N*4)
    
    // Karatsuba 
    Karatsuba32_aux(R, poly1->coefs, poly2->coefs, N1, buf);

    // reduction mod X^N+1
    for (int i = 0; i < N1-1; ++i) result->coefs[i] += R[i] - R[N1+i];
    result->coefs[N1-1] += R[N1-1];
    
    delete[] R;
    delete[] buf;
}


















/* ***********************************************************************************
******************************* Karatsuba 64 *****************************************
*********************************************************************************** */

void torus64PolynomialMultNaive_plain_aux(Torus64* __restrict result, const int* __restrict poly1, const Torus64* __restrict poly2, const int N) {
    const int _2Nm1 = 2*N-1;
    Torus64 ri;
    for (int i=0; i<N; i++) {
        ri=0;
        for (int j=0; j<=i; j++) ri += poly1[j]*poly2[i-j];
        result[i]=ri;
    }
    for (int i=N; i<_2Nm1; i++) {
        ri=0;
        for (int j=i-N+1; j<N; j++) ri += poly1[j]*poly2[i-j];
        result[i]=ri;
    }
}


// A and B of size = size
// R of size = 2*size-1
void Karatsuba64_aux(Torus64* R, const int* A, const Torus64* B, const int size, const char* buf){
    const int h = size / 2;
    const int sm1 = size-1;

    //we stop the karatsuba recursion at h=4, because on my machine,
    //it seems to be optimal
    if (h<=4) {
        torus64PolynomialMultNaive_plain_aux(R, A, B, size);
        return;
    }

    //we split the polynomials in 2
    int* Atemp = (int*) buf; buf += h*sizeof(int);
    Torus64* Btemp = (Torus64*) buf; buf+= h*sizeof(Torus64);
    Torus64* Rtemp = (Torus64*) buf; buf+= size*sizeof(Torus64); 
    //Note: in the above line, I have put size instead of sm1 so that buf remains aligned on a power of 2

    for (int i = 0; i < h; ++i) Atemp[i] = A[i] + A[h+i];
    for (int i = 0; i < h; ++i) Btemp[i] = B[i] + B[h+i];

    // Karatsuba recursivly
    Karatsuba64_aux(R, A, B, h, buf); // (R[0],R[2*h-2]), (A[0],A[h-1]), (B[0],B[h-1])
    Karatsuba64_aux(R+size, A+h, B+h, h, buf); // (R[2*h],R[4*h-2]), (A[h],A[2*h-1]), (B[h],B[2*h-1])
    Karatsuba64_aux(Rtemp, Atemp, Btemp, h, buf);
    R[sm1]=0; //this one needs to be set manually
    for (int i = 0; i < sm1; ++i) Rtemp[i] -= R[i] + R[size+i];
    for (int i = 0; i < sm1; ++i) R[h+i] += Rtemp[i];
}




// poly1, poly2 and result are polynomials mod X^N+1
void torus64PolynomialMultKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomial* poly1, const Torus64Polynomial* poly2, const Globals* env){
    const int N2 = env->N_lvl2;
    Torus64* R = new Torus64[2*N2-1];
    char* buf = new char[32*N2]; //that's large enough to store every tmp variables (2*2*N*8)
    
    // Karatsuba 
    Karatsuba64_aux(R, poly1->coefs, poly2->coefs, N2, buf);

    // reduction mod X^N+1
    for (int i = 0; i < N2-1; ++i) result->coefs[i] = R[i] - R[N2+i];
    result->coefs[N2-1] = R[N2-1];
    
    delete[] R;
    delete[] buf;
}





void torus64PolynomialMultAddKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomial* poly1, const Torus64Polynomial* poly2, const Globals* env){
    const int N2 = env->N_lvl2;
    Torus64* R = new Torus64[2*N2-1];
    char* buf = new char[32*N2]; //that's large enough to store every tmp variables (2*2*N*8)
    
    // Karatsuba 
    Karatsuba64_aux(R, poly1->coefs, poly2->coefs, N2, buf);

    // reduction mod X^N+1
    for (int i = 0; i < N2-1; ++i) result->coefs[i] += R[i] - R[N2+i];
    result->coefs[N2-1] += R[N2-1];
    
    delete[] R;
    delete[] buf;
}




