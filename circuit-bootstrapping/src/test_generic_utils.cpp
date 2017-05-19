#include <cassert>
#include "generic_utils.h"

using namespace std;

class Todo {
    public:
    int a;
    double b;
        Todo(int a, double b): a(a),b(b) {}
};

int main() {
    Todo**** a4 = new_array4<Todo>(10,11,12,13,123456789,7.3);
    for (int i=0; i<10; i++) {
        for (int j=0; j<11; j++) {
            for (int k=0; k<12; k++) {
                for (int l=0; l<13; l++) {
                    assert(a4[i][j][k][l].a==123456789);
                    assert(a4[i][j][k][l].b==7.3);
                    a4[i][j][k][l].a=i+12*j+232*k+90*l;
                }
            }
        }
    }
    for (int i=0; i<10; i++) {
        for (int j=0; j<11; j++) {
            for (int k=0; k<12; k++) {
                for (int l=0; l<13; l++) {
                    assert(a4[i][j][k][l].a==i+12*j+232*k+90*l);
                    assert(a4[i][j][k][l].b==7.3);
                }
            }
        }
    }
    delete_array4(a4);

    ////////////////
    Todo*** a3 = new_array3<Todo>(10,11,12,123456789,7.3);
    for (int i=0; i<10; i++) {
        for (int j=0; j<11; j++) {
            for (int k=0; k<12; k++) {
                assert(a3[i][j][k].a==123456789);
                assert(a3[i][j][k].b==7.3);
                a3[i][j][k].a=i+12*j+232*k;
            }
        }
    }
    for (int i=0; i<10; i++) {
        for (int j=0; j<11; j++) {
            for (int k=0; k<12; k++) {
                assert(a3[i][j][k].a==i+12*j+232*k);
                assert(a3[i][j][k].b==7.3);
            }
        }
    }
    delete_array3(a3);

    ////////////////
    Todo** a2 = new_array2<Todo>(10,11,123456789,7.3);
    for (int i=0; i<10; i++) {
        for (int j=0; j<11; j++) {
                assert(a2[i][j].a==123456789);
                assert(a2[i][j].b==7.3);
                a2[i][j].a=i+12*j;
        }
    }
    for (int i=0; i<10; i++) {
        for (int j=0; j<11; j++) {
                assert(a2[i][j].a==i+12*j);
                assert(a2[i][j].b==7.3);
        }
    }
    delete_array2(a2);

    ////////////////
    Todo* a1 = new_array1<Todo>(10,123456789,7.3);
    for (int i=0; i<10; i++) {
                assert(a1[i].a==123456789);
                assert(a1[i].b==7.3);
                a1[i].a=i;
    }
    for (int i=0; i<10; i++) {
                assert(a1[i].a==i);
                assert(a1[i].b==7.3);
    }
    delete_array1(a1);

}
