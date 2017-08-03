#include <stdint.h>
#include <cmath>

#include <string>
#include <cassert>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <complex>

#include "spqlios-fft.h"

using namespace std;

namespace {

    bool very_close(const complex<double>& a,const complex<double>& b) {
	bool reps = (abs(a-b)<1e-5);
	if (!reps) {
	    cerr << "not close: " << a << " vs. " << b << endl;
	}
	return reps;
    }

    bool very_close(const double& a,const double& b) {
	bool reps = (abs(a-b)<1e-5);
	if (!reps) {
	    cerr << "not close: " << a << " vs. " << b << endl;
	}
	return reps;
    }


    void require(bool condition, const string& message="[condition]") {
	if (!condition) {
	    cerr << "unmet condition: " << condition << endl;
	    abort();
	}
    }
}

int main(int argc, char** argv) {
    static const int nn = 1024;
    void* tables = new_fft_table(nn);
    void* itables = new_ifft_table(nn);
    double* buf_fft = fft_table_get_buffer(tables);
    double* buf_ifft = ifft_table_get_buffer(itables);
    double* a = new double[nn];
    double* a2 = new double[nn];
    double* b = new double[nn];
    for (int i=0; i<nn; i++) a[i]=i;
    for (int i=0; i<nn; i++) a2[i]=i;
    for (int i=0; i<nn; i++) b[i]=a[i];

    for (int i=0; i<nn; i++) buf_ifft[i]=a[i];
    ifft_model(itables);
    for (int i=0; i<nn; i++) a[i]=buf_ifft[i];

    for (int i=0; i<nn; i++) buf_ifft[i]=a2[i];
    ifft(itables,buf_ifft);
    for (int i=0; i<nn; i++) a2[i]=buf_ifft[i];

    for (int i=0; i<nn; i++) 
	require(very_close(a[i],a2[i]));
    
    for (int i=0; i<nn; i++) buf_fft[i]=a[i];
    fft_model(tables);
    for (int i=0; i<nn; i++) a[i]=buf_fft[i];

    for (int i=0; i<nn; i++) buf_fft[i]=a2[i];
    fft(tables,buf_fft);
    for (int i=0; i<nn; i++) a2[i]=buf_fft[i];

    for (int i=0; i<nn; i++) 
	require(very_close(a[i],a2[i]));
    for (int i=0; i<nn; i++) 
	require(very_close(a[i],nn/2*b[i]));

#ifdef NDEBUG
    int nbtrials = 100000;
    long begin = clock();
    for (int trials=0; trials<nbtrials; trials++)
	fft(tables,buf_fft);
    long end = clock();
    cout << "fft_time " << double(end-begin)/nbtrials << endl;
    begin = clock();
    for (int trials=0; trials<nbtrials; trials++)
	ifft(tables,buf_ifft);
    end = clock();
    cout << "ifft_time " << double(end-begin)/nbtrials << endl;
    begin = clock();
    for (int trials=0; trials<nbtrials; trials++)
	fft_model(tables);
    end = clock();
    cout << "fft_model_time " << double(end-begin)/nbtrials << endl;
    begin = clock();
    for (int trials=0; trials<nbtrials; trials++)
	ifft_model(tables);
    end = clock();
    cout << "ifft_model_time " << double(end-begin)/nbtrials << endl;
#endif
}


