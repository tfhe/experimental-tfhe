#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <thread>
#include <vector>
#include <omp.h>

using namespace std;

void f(int j) {
    static __thread int threadnumber;
    string s = "Hello World;";
#pragma omp critical
    {
        threadnumber = j;
        for (int i = 0; i < s.size(); i++) {
            cerr << s[i];
            usleep(10000);
        }
        cerr << "begin " << threadnumber << endl;
    }
    sleep(1);
#pragma omp critical
    {
        for (int i = 0; i < s.size(); i++) {
            cerr << s[i];
            usleep(10000);
        }
        cerr << "end " << threadnumber << endl;
    }
}


int main1() {
    vector < thread * > toto(16);

    for (int i = 0; i < 16; i++) {
        toto[i] = new thread(f, i);
    }

    for (int i = 0; i < 16; i++) {
        toto[i]->join();
    }


}


int main() {

#pragma omp parallel for
    for (int i = 0; i < 16; i++) {
        f(i);
    }

}
