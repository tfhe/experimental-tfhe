#ifndef GENERIC_UTILS_H
#define GENERIC_UTILS_H

#include <random>

//new array and multidimensional arrays
template<typename T, typename... U>
T* new_array1(int a, U... params) {
    size_t s = 1*sizeof(int) + a*sizeof(T);
    char* chunk = (char*) malloc(s);
    int* sp = (int*) chunk; 
    T* data = (T*) (chunk + 1*sizeof(int));
    sp[0] = a;
    for (int i=0; i<a; ++i)
        new (data+i) T(params...);
    return data;
};

template<typename T>
void delete_array1(T* data) {
    char* chunk = (char*) data;
    int* sp = (int*) (chunk - 1*sizeof(int));
    int a = sp[0];
    for (int i=0; i<a; ++i)
        (data+i)->~T();
    free(sp); //address of the original chunk
};


template<typename T, typename... U>
T** new_array2(int a, int b, U... params) {
    size_t s = 2*sizeof(int) + a*b*sizeof(T) + a*sizeof(T*);
    char* chunk = (char*) malloc(s);
    int* sp = (int*) chunk; 
    T** data2 = (T**) (chunk + 2*sizeof(int));
    T* data1 = (T*) (chunk + 2*sizeof(int) + a*sizeof(T*)); 
    sp[0] = a;
    sp[1] = b;
    for (int i=0; i<a*b; ++i)
        new (data1+i) T(params...);
    for (int i=0; i<a; ++i)
        data2[i]=data1+i*b;
    return data2;
};

template<typename T>
void delete_array2(T** data) {
    char* chunk = (char*) data;
    int* sp = (int*) (chunk - 2*sizeof(int));
    int a = sp[0];
    int b = sp[1];
    T* data1 = data[0]; 
    for (int i=0; i<a*b; ++i)
        (data1+i)->~T();
    free(sp); //address of the original chunk
};


template<typename T, typename... U>
T*** new_array3(int a, int b, int c, U... params) {
    size_t s = 3*sizeof(int) + a*b*c*sizeof(T) + a*b*sizeof(T*) + a*sizeof(T**);
    char* chunk = (char*) malloc(s);
    int* sp = (int*) chunk; 
    T*** data3 = (T***) (chunk + 3*sizeof(int));
    T** data2 = (T**) (chunk + 3*sizeof(int) + a*sizeof(T**)); 
    T* data1 =  (T*) (chunk + 3*sizeof(int) + a*sizeof(T**) + a*b*sizeof(T*)); 
    sp[0] = a;
    sp[1] = b;
    sp[2] = c;
    for (int i=0; i<a*b*c; ++i)
        new (data1+i) T(params...);
    for (int i=0; i<a*b; ++i)
        data2[i]=data1+i*c;
    for (int i=0; i<a; ++i)
        data3[i]=data2+i*b;
    return data3;
};

template<typename T>
void delete_array3(T*** data) {
    char* chunk = (char*) data;
    int* sp = (int*) (chunk - 3*sizeof(int));
    int a = sp[0];
    int b = sp[1];
    int c = sp[2];
    T* data1 = data[0][0]; 
    for (int i=0; i<a*b*c; ++i)
        (data1+i)->~T();
    free(sp); //address of the original chunk
};

template<typename T, typename... U>
T**** new_array4(int a, int b, int c, int d, U... params) {
    size_t s = 4*sizeof(int) + a*b*c*d*sizeof(T) + a*b*c*sizeof(T*) + a*b*sizeof(T**) + a*sizeof(T***);
    char* chunk = (char*) malloc(s);
    int* sp = (int*) chunk; 
    T**** data4 = (T****) (chunk + 4*sizeof(int));
    T*** data3 = (T***) (chunk + 4*sizeof(int) + a*sizeof(T***));
    T** data2 = (T**) (chunk + 4*sizeof(int) + a*sizeof(T***) + a*b*sizeof(T**)); 
    T* data1 = (T*) (chunk + 4*sizeof(int) + a*sizeof(T***) + a*b*sizeof(T**)+ a*b*c*sizeof(T*)); 
    sp[0] = a;
    sp[1] = b;
    sp[2] = c;
    sp[3] = d;
    for (int i=0; i<a*b*c*d; ++i)
        new (data1+i) T(params...);
    for (int i=0; i<a*b*c; ++i)
        data2[i]=data1+i*d;
    for (int i=0; i<a*b; ++i)
        data3[i]=data2+i*c;
    for (int i=0; i<a; ++i)
        data4[i]=data3+i*b;
    return data4;
};

template<typename T>
void delete_array4(T**** data) {
    char* chunk = (char*) data;
    int* sp = (int*) (chunk - 4*sizeof(int));
    int a = sp[0];
    int b = sp[1];
    int c = sp[2];
    int d = sp[3];
    T* data1 = data[0][0][0]; 
    for (int i=0; i<a*b*c*d; ++i)
        (data1+i)->~T();
    free(sp); //address of the original chunk
};

//////////////
//#define FALSE_RANDOM
#ifdef FALSE_RANDOM

inline int random_bit() { return 1; }
inline int32_t random_int32() { return 0xcccccccc; }
inline int64_t random_int64() { return 0xccccccccccccccccul; }
inline double random_gaussian_double(double center, double stdev) {
    return center; 
}

inline int32_t random_gaussian32(int32_t center, double stdev) {
    return center;    
}

inline int64_t random_gaussian64(int64_t center, double stdev) {
    return center;    
}

class Random {};

#else
/////// TRUE RANDOM
class Random {
    public:
        std::default_random_engine generator;
        std::uniform_int_distribution<int> bit_distribution;
        std::uniform_int_distribution<int32_t> int32_distribution;
        std::uniform_int_distribution<int64_t> int64_distribution;
        std::normal_distribution<double> gaussian_distribution;
        Random():
            bit_distribution(0,1),
            int32_distribution(std::numeric_limits<int32_t>::min(),std::numeric_limits<int32_t>::max()),
            int64_distribution(std::numeric_limits<int64_t>::min(),std::numeric_limits<int64_t>::max()),
            gaussian_distribution(0.,1.) {}
};

extern Random* global_random; 

inline int random_bit() { return global_random->bit_distribution(global_random->generator); }
inline int32_t random_int32() { return global_random->int32_distribution(global_random->generator); }
inline int64_t random_int64() { return global_random->int64_distribution(global_random->generator); }
inline double random_gaussian_double(double center, double stdev) { 
    return stdev*global_random->gaussian_distribution(global_random->generator)+center; 
}

inline int32_t random_gaussian32(int32_t center, double stdev) {
    static const double _2p32 = pow(2.,32);
    double val = stdev*global_random->gaussian_distribution(global_random->generator)*_2p32;
    int32_t ival = (int32_t) val;
    return ival+center;    
}

inline int64_t random_gaussian64(int64_t center, double stdev) {
    static const double _2p64 = pow(2.,64);
    double val = stdev*global_random->gaussian_distribution(global_random->generator)*_2p64;
    int64_t ival = (int64_t) val;
    //printf("ival: %ld\n", ival);
    return ival+center;    
}
#endif
////////////

#endif // GENERIC_UTILS_H
