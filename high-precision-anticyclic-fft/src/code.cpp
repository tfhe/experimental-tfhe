#include <complex>
#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>
#include <random>

using namespace std;

typedef int64_t Torus64;

typedef int32_t Torus32;

const double DBL2M64 = 1./(1l<<32)/(1l<<32);

#if 1
//fixed point real number:
//96 bit integer v
//representing the value v/2^64 mod 2^32

#ifndef intmul
#define intmul intmul_best
#endif

struct Real96 {
    __uint128_t v;

    Real96(uint64_t value=0);
    Real96(uint64_t lo, uint64_t hi);
};

Real96::Real96(uint64_t value) {
    v=value;
}

Real96::Real96(uint64_t lo, uint64_t hi) {
    v=hi;
    v<<=64;
    v|=lo;
}

ostream& operator<<(ostream& out, const Real96& r) {
    out << "["<< (int64_t(r.v>>64)) << "].[" << (int64_t(r.v&0xFFFFFFFFFFFFFFFFUL)) << "] (" << double(__int128_t(r.v))*DBL2M64 << ")";
    return out;
}

void add(Real96& dest, const Real96& a, const Real96& b) {
    dest.v=a.v+b.v;
}
Real96 operator+(const Real96& a, const Real96& b) {
    Real96 reps; add(reps,a,b); return reps;
}
void operator+=(Real96& a, const Real96& b) {
    add(a,a,b);
}


void sub(Real96& dest, const Real96& a, const Real96& b) {
    dest.v=a.v-b.v;
}
void neg(Real96& dest, const Real96& a) {
    dest.v=-a.v;
}
Real96 operator-(const Real96& a, const Real96& b) {
    Real96 reps; sub(reps,a,b); return reps;
}
void operator-=(Real96& a, const Real96& b) {
    sub(a,a,b);
}
Real96 operator-(const Real96& a) {
    Real96 reps; neg(reps,a); return reps;
}

void extmul(Real96& dest, int a, const Real96& b) {
    dest.v=a*b.v;
}


void  intmul_ref(Real96& dest, const Real96& a, const Real96& b) {
    __int128_t av = a.v;
    //__uint128_t avu = a.v;
    __int128_t bv = b.v;
    //__uint128_t bvu = b.v;
    //__uint128_t alo = av&0xFFFFFFFFFFFFFFFFUL;
    __uint128_t blo = bv&0xFFFFFFFFFFFFFFFFUL;
    int32_t ahi = av>>64;
    __int128_t bhi = bv>>127;
    assert(bhi==0 || bhi==-1);
    assert(ahi<1l<<31 || ahi>=-1l<<31);

    __int128_t w0 = (uint64_t(a.v)*blo)>>64; //between [0 and 2^64-2]
    w0 += ahi*blo - (bhi&av);
    dest.v = w0;

}

void __attribute ((noinline)) intmul_asm_fk(Real96& dest, const Real96& a, const Real96& b) {
    //const Real96* aad = &a;
    //const Real96* bad = &b;
    //Real96* destad = &dest;
    uint64_t tmpa;
    uint64_t tmpb;
    uint64_t tmpc;
    uint64_t tmpd;
    uint64_t tmpe;
    uint64_t tbval;
    __asm volatile (
	"movq      %[b],%[tbval]\n"
        "movq      8(%[a]), %[tmpb]\n"           // tmpb=ahi
        "sarq      $63, %[tmpb]\n"               // tmpb=ahi<0?-1:0
        "movq      (%[tbval]), %%rdx\n"          // rdx=blo
        "mulx      (%[a]), %[tmpa], %[tmpe]\n"   // [e:a] = alo*blo = [u:v]
        "mulx      8(%[a]), %[tmpa], %[tmpd]\n"  // [d:a] = ahi*blo = [w:x]
        "negq      %%rdx\n"
        "andq      %%rdx, %[tmpb]\n"             // tmpb=ahi<0?-blo:0
        "addq      %[tmpe], %[tmpa]\n"           // x+=u
        "movq      8(%[tbval]), %[tmpc]\n"       // c=bhi
        "adcq      %[tmpb], %[tmpd]\n"           // w+=-blo if ahi<0
        "btq       $63, %[tmpc]\n"               // test sign of bhi
        "jnc       1f\n"
        "subq      (%[a]), %[tmpa]\n"            // x-=alo
        "sbbq      8(%[a]), %[tmpd]\n"           // w-=ahi+carry
        "1:\n"
        "movq      %[tmpa], (%[dest])\n"
        "movq      %[tmpd], 8(%[dest])\n" 
	: [tmpa] "=r" (tmpa),
	  [tmpb] "=r" (tmpb),
	  [tmpc] "=r" (tmpc),
	  [tmpd] "=r" (tmpd),
	  [tmpe] "=r" (tmpe),
	  [tbval] "=r" (tbval)
	: [dest] "r" (&dest),
	  [a] "r"(&a),
	  [b] "r"(&b)
	: "%rdx","memory"
	   );
}

#include <bmi2intrin.h>

typedef union {
    __uint128_t vu128;
    unsigned long long vu64[2];
    __int128_t v128;
    long long v64[2];
} UINT128;

void intmul_best(Real96& dest, const Real96& a, const Real96& b) {
    const UINT128* aa = (const UINT128*) &a;
    const UINT128* ab = (const UINT128*) &b;
    //UINT128* adest = (UINT128*) &dest;
    UINT128 tab;
    UINT128 tcd;
    //UINT128 tef;
    _mulx_u64(ab->vu64[0],aa->vu64[0],&tcd.vu64[0]);  // td=u
    tab.vu64[0]=_mulx_u64(ab->vu64[0],aa->vu64[1],&tab.vu64[1]);
    tcd.vu64[1]=(aa->v64[1]>>63)&(-ab->vu64[0]); //(aa->vu64[1]&0x8000000000000000UL)?(-ab->vu64[0]):0;
    tab.vu128+=tcd.vu128;
    if (ab->vu64[1]&0x8000000000000000UL) tab.vu128-=aa->vu128;
    //tcd.vu64[0]=ab->vu64[1]&aa->vu64[0];
    //tcd.vu64[1]=ab->vu64[1]&aa->vu64[1];
    //tab.vu128-=tcd.vu128;
#ifdef NDEBUG
    dest.v=tab.vu128;
#else
    intmul_ref(dest,a,b);
    assert(dest.v==tab.vu128);
#endif
}

Real96 operator*(int a, const Real96& b) {
    Real96 reps; extmul(reps,a,b); return reps;
}
void operator*=(Real96& b, int a) {
    extmul(b,a,b);
}
Real96 operator*(const Real96& a, const Real96& b) {
    Real96 reps; intmul(reps,a,b); return reps;
}
void operator*=(Real96& a, const Real96& b) {
    intmul(a,a,b);
}

Real96 t64tor96(Torus64 v) {
    Real96 reps; 
    int64_t vv = v;
    reps.v = vv;
    return reps;
}

Real96 dtor96(double v) {
    Real96 reps;
    int64_t iv = floor(v);
    v = v-iv;
    assert(v>=0 && v<1.);
    v *= (1ul<<32);
    v *= (1ul<<32);
    uint64_t lv = v;
    __int128_t pv = iv;
    pv <<= 64;
    pv |= lv;
    reps.v = pv;
    return reps;
}

bool very_close(const Real96& a,const Real96& b) {
    bool reps = (abs(__int128_t(a.v-b.v))<10000);
    if (!reps) {
	cerr << "not close: " << a << " vs. " << b << endl;
    }
    return reps;
}

#else
typedef double Real96;

Real96 t64tor96(Torus64 v) {
    return double(v)*DBL2M64;
}

Real96 dtor96(double v) {
    return v;
}

bool very_close(const Real96& a,const Real96& b) {
    bool reps = (abs(a-b)<1e-5);
    if (!reps) {
	cerr << "not close: " << a << " vs. " << b << endl;
    }
    return reps;
}

#endif

//-----------------------------------------------------------------------

typedef complex<Real96> Cplx96;

//-----------------------------------------------------------------------

#include <NTL/RR.h>
using namespace NTL;

RR PI_RR = ComputePi_RR();

Real96 accurate_cos(int i,int n) { //cos(2pi*i/n)
    i = ((i%n) + n)%n;
    if (i==0) return Real96(-1,0);
    ZZ cosi = RoundToZZ(cos(PI_RR*2*i/n)*pow(to_RR(2),to_RR(64)));
    if (cosi>=0)
	return Real96(to_long(cosi));
    else
	return Real96(to_long(cosi+power(to_ZZ(2),64)),-1);

/*
    if (i>=3*n/4) return dtor96(cos(2.*M_PI*(n-i)/double(n)));
    if (i>=2*n/4) return dtor96(-cos(2.*M_PI*(i-n/2)/double(n)));
    if (i>=1*n/4) return dtor96(-cos(2.*M_PI*(n/2-i)/double(n)));
    return dtor96(cos(2.*M_PI*(i)/double(n)));
    */
}

Real96 accurate_sin(int i,int n) { //sin(2pi*i/n)
    i = ((i%n) + n)%n;
    if (i==n/4) return Real96(-1,0);
    ZZ sini = RoundToZZ(sin(PI_RR*2*i/n)*pow(to_RR(2),to_RR(64)));
    if (sini>=0)
	return Real96(to_long(sini));
    else
	return Real96(to_long(sini+power(to_ZZ(2),64)),-1);
/*
    if (i>=3*n/4) return dtor96(-sin(2.*M_PI*(n-i)/double(n)));
    if (i>=2*n/4) return dtor96(-sin(2.*M_PI*(i-n/2)/double(n)));
    if (i>=1*n/4) return dtor96(sin(2.*M_PI*(n/2-i)/double(n)));
    return dtor96(sin(2.*M_PI*(i)/double(n)));
    */
}

//reverse the bits of i (mod n)
int rev(int i, int n) {
    int reps=0;
    for (int j=1; j<n; j*=2) {
	reps = 2*reps + (i%2);
	i/=2;
    }
    return reps;
}

bool very_close(const Cplx96& a,const Cplx96& b) {
    bool reps = (very_close(a.real(),b.real()) && very_close(a.imag(),b.imag()));
    if (!reps) {
	cerr << "not close: " << a << " vs. " << b << endl;
    }
    return reps;
}






// FFT from Torus64^N to Cplx96^(N/2)  mod X^N+1
// N = 2048 (note: n=2N ici)


//at the beginning of iteration nn
// a_{j,i} has P_{i%nn}(omega^j) 
// where j between [rev(1) and rev(3)[
// and i between [0 and nn[
void ifft_check(int n, int nn, const Cplx96* acur, const vector<Real96>& a, const vector<Cplx96>& powomega) {
    int ns4=n/4;	
    cerr << "Checking iteration " << nn << endl;
    for (int i=0; i<ns4; i++) {
	cout << "i: " << i << "   " << acur[i] << endl;
    }
    int m = n/nn;
    int rev1m = rev(1,m);
    int rev3m = rev(3,m);
    int idex = 0;
    for (int revj=rev1m; revj<rev3m; revj++) {
	int j = rev(revj,m);
	cerr << "check-- j: " << j << endl;
	for (int i=0; i<nn; i++) {
	    cerr << "check--- i: " << i << "(mod " << nn << ")" << endl;
	    const Cplx96& test_cur = acur[idex];
	    //sum_[t=i%nn] a_t omega^jt
	    Cplx96 pij(0,0);
	    for (int k=i; k<n; k+=nn) {
		//cout << "ak: " << a[k] << endl;
		pij += a[k] * powomega[(k*j) % n];
	    }
	    assert(very_close(test_cur,pij));
	    idex++;
	}
    }
}


//at the beginning of iteration halfnn:
//   m=n/halfnn
//   P_{j%m}(omb^i)
//   for j in [rev(1,m) to rev(3,m)[
//   and i in [0,halfnn[
void fft_check(
	int n, int halfnn, 
	const Cplx96* pcur,
	const vector<Cplx96>& p,
	const vector<Cplx96>& powombar
	) {
    int ns4=n/4;
    cerr << "DIRECT FFT: Checking iteration " << halfnn  << endl;
    for (int i=0; i<ns4; i++) {
	cout << "i: " << i << "   " << pcur[i] << endl;
    }
    int m = n/halfnn;
    int rev1m = rev(1,m);
    int rev3m = rev(3,m);
    int idex = 0;
    for (int revj=rev1m; revj<rev3m; revj++) {
	int j = rev(revj,m);
	cerr << "check-- j: " << j << "(mod " << m << ")" << endl;
	for (int i=0; i<halfnn; i++) {
	    cerr << "check--- i: " << i << endl;
	    //P_sum_[k=j%m] p_k omb^ik-j
	    Cplx96 pij(0,0);	
	    for (int k=j; k<n; k+=m) {
		//if (halfnn==8 && j==1 && i==1) cerr << "pij(" << pij << ")" << "+= p_"<<k<<"("<<p[k]<<") * omb["<<i*(k-j)<<"]("<< powombar[(i*(k-j)) % n] <<")" << endl; 
		pij += p[k] * powombar[(i*(k-j)) % n];
	    }
	    assert(very_close(pcur[idex],pij));
	    idex++;
	}
    }
}



void precomp_iFFT(vector<Cplx96>& powomega, int n) {
    powomega.resize(n);
    for (int i=0; i<n; i++)
	powomega[i]=Cplx96(accurate_cos(i,n),accurate_sin(i,n));
}

void precomp_FFT(vector<Cplx96>& powombar, int n) {
    powombar.resize(n);
    for (int i=0; i<n; i++)
	powombar[i]=Cplx96(accurate_cos(i,n),accurate_sin((n-i)%n,n));
}

// P -> P(omega)
void iFFT(Cplx96* out, const Torus64* in, int n, const vector<Cplx96>& powomega) {
    //const int N = n/2;
    const int ns4 = n/4;

#ifndef NDEBUG
    vector<Real96> a; a.resize(n);
    for (int i=0; i<n/2; i++)
	a[i]=t64tor96(in[i]/2);
    for (int i=0; i<n/2; i++)
	a[n/2+i]=-a[i];
#endif


    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j=0; j<ns4; j++)
	out[j] = Cplx96(t64tor96(in[j]),t64tor96(in[j+ns4]))*powomega[j];

    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j) 
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn=ns4; nn>=2; nn/=2) {
	int halfnn = nn/2;
#ifndef NDEBUG
	cerr << "Starting iteration " << nn << endl;
	int m = n/nn;
	ifft_check(n, nn, out, a, powomega);
#endif
	for (int block=0; block<ns4; block+=nn) {
#ifndef NDEBUG
	    int j = rev(rev(1,m)+block,m);
	    cerr << "-- block j: " << j  << " --> " << j << "," << j+halfnn/2 << endl;
#endif
	    for (int off=0; off<halfnn; off++) {
#ifndef NDEBUG
		cerr << "--- i: " << off << " using: omg^" << (2*(ns4/halfnn)*off)%n << endl;
#endif
		Cplx96 t1 = out[block+off];
		Cplx96 t2 = out[block+off+halfnn];
		out[block+off]=t1+t2;
		out[block+off+halfnn]=(t1-t2)*powomega[(2*(ns4/halfnn)*off)%n];
	    }
	}
    }
    {
#ifndef NDEBUG
	int nn = 1;
	ifft_check(n, nn, out, a, powomega);
#endif
    }
}

// P(omega) -> P
void FFT(Torus64* out, Cplx96* in, int n, const vector<Cplx96>& powombar) {
    //const int N = n/2;
    const int ns4 = n/4;

#ifndef NDEBUG
    vector<Cplx96> a; a.resize(n);
    for (int i=0; i<n; i++) a[i]=Cplx96(0);
    int rev1m = rev(1,n);
    int rev3m = rev(3,n);
    for (int revj=rev1m; revj<rev3m; revj++) {
	int j = rev(revj,n);
	cout << "assign:" << j << " " << revj-rev1m << endl;
	a[j]=in[revj-rev1m];
	a[n-j]=Cplx96(a[j].real(),-a[j].imag());
    }
#endif


#ifndef NDEBUG
	cerr << "Checking iteration 1" << endl;
	fft_check(n, 1, in, a, powombar);
#endif

    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j) 
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn=2; nn<=ns4; nn*=2) {
	int halfnn = nn/2;
	for (int block=0; block<ns4; block+=nn) {
	    //#ifndef NDEBUG
	    //	    int j = rev(rev(1,m)+block,m);
	    //	    cerr << "-- block j: " << j  << " --> " << j << "," << j+halfnn/2 << endl;
	    //#endif
	    for (int off=0; off<halfnn; off++) {
		//#ifndef NDEBUG
		//		cerr << "--- i: " << off << " using: omg^" << (2*(ns4/halfnn)*off)%n << endl;
		//#endif
		Cplx96 t1 = in[block+off];
		Cplx96 t2 = in[block+off+halfnn]*powombar[(2*(ns4/halfnn)*off)%n];
		in[block+off]=t1+t2;
		in[block+off+halfnn]=(t1-t2);
	    }
	}
#ifndef NDEBUG
	cerr << "Ending iteration " << nn << endl;
	fft_check(n, nn, in, a, powombar);
#endif
    }


    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j=0; j<ns4; j++) {
	in[j] *= powombar[j];
	out[j] = in[j].real().v>>10;  // /ns4;  //divide by N/2
	out[j+ns4] = in[j].imag().v>>10; // /ns4; //divide by N/2
    }

    {
	//#ifndef NDEBUG
	//	int nn = 1;
	//	ifft_check(n, nn, out, a, powomega);
	//#endif
    }
}



int main(int argc, char** argv) {

#if 1

    const int N = 2048;
    const int Ns2 = N/2;
#ifndef NDEBUG
    const int NBTRIALS=1;
#else
    const int NBTRIALS=10000;
#endif

#ifndef NDEBUG
    //test cosinus
    for (int i=0; i<=2*N; i++) {
	//if (i!=1024) continue;
	Real96 c = accurate_cos(i,2*N);
	Real96 s = accurate_sin(i,2*N);
	cout << "i: " << i  << ", " << cos(2*i*M_PI/64.) << endl;
	cout << "c " << c << endl;
	cout << "s " << s << endl;
	cout << "c*c " << c*c << endl;
	cout << "s*s " << s*s << endl;
	Real96 t = c*c+s*s;
	cout << "cos test" <<  i << ":" << t << endl;
	assert(very_close(t,Real96(0,1)));
    }
#endif

    //test de la IFFT
    Cplx96 out[Ns2];
    Torus64 in[N];
    Torus64 revout[N];

    vector<Cplx96> powomega;
    vector<Cplx96> powombar;

    std::default_random_engine generator;
    std::uniform_int_distribution<Torus64> distribution(numeric_limits<Torus64>::min(),numeric_limits<Torus64>::max());

    for (int i=0; i<N; i++) {
	in[i]=distribution(generator);
    }

    precomp_iFFT(powomega,2*N);
    precomp_FFT(powombar,2*N);

#ifndef NDEBUG
    //test powomega
    for (int i=0; i<2*N; i++) {
	//if (i!=1024) continue;
	cout << "powomega: " << i << " : " << powomega[i] << endl; 
	assert(very_close(powomega[i]*powomega[(2*N-i)%(2*N)],Cplx96(Real96(0,1))));
	assert(very_close(powomega[i]*powombar[i],Cplx96(Real96(0,1))));
    }

#endif

    clock_t t0 = clock();
    for (int i=0; i<NBTRIALS; i++)
	iFFT(out,in,2*N,powomega);
    clock_t t1 = clock();
    for (int i=0; i<NBTRIALS; i++)
	FFT(revout,out,2*N,powombar);
    clock_t t2 = clock();
#ifndef NDEBUG
    for (int i=0; i<N; i++)
	cout << hex << revout[i] << " "<< in[i] << " " << (revout[i]^in[i]) << endl;
#endif
    cout << "time IFFT: " << (t1-t0)/double(NBTRIALS) << "mus" << endl;
    cout << "time FFT: " << (t2-t1)/double(NBTRIALS) << "mus" << endl;

#endif

#if 0
    const int NBTRIALS=100;
    
    //test d'un keyswitch binaire
    int nlvl2=2048;
    int tlvl2=32;
    int klvl1=2;
    int Nlvl1=1024;
    int llvl1=4;

    int KS_size=Nlvl1*(klvl1+1)*(klvl1+1)*tlvl2*nlvl2;
    cout << "KS_size: " << KS_size << " * 32bits" << endl;
    Torus32* KS_rrr = new Torus32[Nlvl1*(klvl1+1)*(klvl1+1)*tlvl2*nlvl2];
    Torus32** KS_rr = new Torus32*[(klvl1+1)*tlvl2*nlvl2];
    Torus32*** KS_r = new Torus32**[tlvl2*nlvl2];
    Torus32**** KS = new Torus32***[nlvl2];
    for (int i=0; i<Nlvl1*(klvl1+1)*(klvl1+1)*tlvl2*nlvl2; i++) KS_rrr[i]=i; //random
    for (int i=0; i<(klvl1+1)*tlvl2*nlvl2; i++) KS_rr[i]=KS_rrr+i*(Nlvl1*(klvl1+1));
    for (int i=0; i<tlvl2*nlvl2; i++) KS_r[i]=KS_rr+i*(klvl1+1);
    for (int i=0; i<nlvl2; i++) KS[i]=KS_r+i*tlvl2;

    Torus64* samples_r = new Torus64[(nlvl2+1)*llvl1];
    Torus64** samples = new Torus64*[llvl1];
    for (int i=0; i<(nlvl2+1)*llvl1; i++) samples_r[i]=i; //random
    for (int i=0; i<llvl1; i++) samples[i]=samples_r+i*(nlvl2+1);

    Torus32* KSres = new Torus32[(klvl1+1)*Nlvl1];    

    clock_t tdebks = clock();
    for (int trial=0; trial<NBTRIALS; trial++) {
	int u=0;
	for (int i=0; i<(klvl1+1)*Nlvl1; i++) KSres[i]=0;
	KSres[u*Nlvl1]=samples[0][nlvl2]; //b
	for (int i=0; i<nlvl2; i++) {
	    //decompose a
	    uint64_t ai=samples[0][i];
	    for (int j=0; j<tlvl2; j++) {
		int8_t aij=(ai & 1ul<<(63-j))?1:0;
		for (int l=0; l<(klvl1+1)*Nlvl1; l++)
		    KSres[l]+=aij*KS[i][j][u][l];
	    }
	}
    }
    clock_t tendks = clock();

    cout << "KS:" << double(tendks-tdebks)/NBTRIALS << "mus " << endl;
#endif

}
