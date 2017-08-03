#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <thread>
#include <vector>
#include <omp.h>
#include <poc_types.h>

using namespace std;


void torus32PolynomialMultNaive_plain_aux(Torus32* __restrict result, const int* __restrict poly1, const Torus32* __restrict poly2, const int N);
void torus32PolynomialMultKaratsuba_lvl1(Torus32Polynomial* result, const IntPolynomial* poly1, const Torus32Polynomial* poly2, const int N1);
// void torusPolynomialMultFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);





// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}






int main(int argc, char **argv) {
    const int nb_products = 10; // number of times we test the products
    const int N = 4;

    // input polynomials
    IntPolynomial* a = new IntPolynomial(N);
    int* a_coefs = new int[N];
    Torus32Polynomial* b = new Torus32Polynomial(N);
    Torus32* b_coefs = new Torus32[N];
    //TorusPolynomial* b1 = new TorusPolynomial(N);

    // output polynomials
    Torus32* resNaive_coefs = new Torus32[N];
    Torus32Polynomial* resKaratsuba = new Torus32Polynomial(N);
    //TorusPolynomial* resFFT = new TorusPolynomial(N);


    // counting the clock cycles
    double nb_cycles_naive[nb_products];
    double nb_cycles_karatsuba[nb_products];
    //double nb_cycles_fft[nb_products];


    // Parallelization of the multiplications with openmp
//#pragma omp parallel for
    for (int i = 0; i < nb_products; ++i) {
        cout << "trial " << i << endl;
        for (int j = 0; j < N; j++) {
            a->coefs[j] = (rand() % (4 * N)) - (2 * N);
            a_coefs[j] = a->coefs[j];
            b->coefs[j] = rand();
            b_coefs[j] = b->coefs[j];
            //b1->coefsT[j] = b->coefs[j];

            resNaive_coefs[j] = 0;
            resKaratsuba->coefs[j] = 0;
            //resFFT->coefsT[j] = 0;
        }
        cout << "Input polynomials: generated" << endl;
        cout << "Output polynomials: initialized to 0" << endl;

        //measure the execution time
        clock_t cstart, cend;

        // Naive multiplication
        cstart = clock();
        //torusPolynomialMultNaive(resNaive, a, b);
        torus32PolynomialMultNaive_plain_aux(resNaive_coefs, a_coefs, b_coefs, N);
        /*
        for (int j = 0; j < N; ++j) {
            cout << resNaive_coefs[j] << " " ;
        }
        cout << endl;
        */
        cout << "Naive multiplication: done" << endl;
        cend = clock();
        nb_cycles_naive[i] = cend - cstart;


        // Karatsuba multiplication
        cstart = clock();
        //torusPolynomialMultKaratsuba(resKaratsuba, a, b);
        torus32PolynomialMultKaratsuba_lvl1(resKaratsuba, a, b, N);
        /*
        for (int j = 0; j < N; ++j) {
            cout << resKaratsuba->coefs[j] << " " ;
        }
        cout << endl;
        */
        cout << "Karatsuba multiplication: done" << endl;
        cend = clock();
        nb_cycles_karatsuba[i] = cend - cstart;


        /*
        // FFT multiplication
        cstart = clock();
        torusPolynomialMultFFT(resFFT, a, b1);
        cout << "FFT multiplication: done" << endl;
        //TorusPolynomial_ifft(test_fft,resNaive);
        //TorusPolynomial_fft(resFFT,test_fft);
        cend = clock();
        nb_cycles_fft[i] = cend - cstart;
        */


        for (int j = 0; j < N; j++) {
            if (abs(resNaive_coefs[j] - resKaratsuba->coefs[j]) > 1) {
                cout << abs(resNaive_coefs[j] - resKaratsuba->coefs[j]) << endl;
                cerr << j << " " << resNaive_coefs[j] << " vs. " << resKaratsuba->coefs[j] << endl;
                //dieDramatically("Naive != Karatsuba\n");
            }
        }


        /*
       for (int j = 0; j < N; j++) {
            if (abs(int(resNaive_coefs[j] - resFFT->coefs[j])) > 1) {
                cerr << j << " " << resNaive_coefs[j] << " vs. " << resFFT->coefs[j] << endl;
                dieDramatically("Naive != FFT\n");
        }
        */

    }

    // computing the average number of cycles per type of multiplication
    double cyc_naive = 0;
    double cyc_karatsuba = 0;
    //double cyc_fft = 0;
    for (int i = 0; i < nb_products; ++i) {
        cyc_naive += nb_cycles_naive[i];
        cyc_karatsuba += nb_cycles_karatsuba[i];
        //cyc_fft += nb_cycles_fft[i];
    }
    cyc_naive = cyc_naive / nb_products;
    cyc_karatsuba = cyc_karatsuba / nb_products;
    //cyc_fft = cyc_fft / nb_products;

    cout << "torusPolynomialMultNaive: " << cyc_naive << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultNaive time: " << (cyc_naive / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    cout << "torusPolynomialMultKaratsuba: " << cyc_karatsuba << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultKaratsuba time: " << (cyc_karatsuba / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    /*
    cout << "torusPolynomialMultFFT: " << cyc_fft << " clock cycles (average)" << endl;
    cout << "torusPolynomialMultFFT time: " << (cyc_fft / CLOCKS_PER_SEC) << " seconds" << endl;
    cout << endl;
    */


/*
    delete_IntPolynomial(a);
    delete_TorusPolynomial(b);
    delete_TorusPolynomial(resNaive);
    delete_TorusPolynomial(resKaratsuba);
    //delete_TorusPolynomial(resFFT);
*/
    return 0;
}
