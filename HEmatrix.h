#pragma once
#include "common.h"

class HEvec {
public:
    Context context;
    long logslot;
    long slot;
    Ciphertext* ct;
    long dim;

    HEvec(Context &context) {
        this->context = context;
        this->logslot = getLogFullSlots(context);
        this->slot = 1 << logslot;
    }

    void setHEvec(Ciphertext* &ct, long dim) {
        this->dim = dim;
        if(dim < slot) {
            this->ct = ct;
        }
        else {
            cout << "아직 구현 안됬습니다 ^__^\n";
        }

    }
};

// Assumtion : any HEmatrix is square matrix
class HEmatrix {
public: 
    Context context;
    long logslot;
    long slot;
    vector<Ciphertext*> mat;
    long matsize;

    HEmatrix(Context &context) {
        this->context = context;
        this->logslot = getLogFullSlots(context);
        this->slot = 1 << logslot;
    };

    void setHEmatrix(vector<Ciphertext*> &ctxts, long size) {
        matsize = size;
        if(size < slot) {
            mat.resize(size);
            for(int i=0; i<size; ++i) {
                mat[i] = new Ciphertext(context);
                mat[i] = ctxts[i];
            }
        }
        else {
            cout << "아직 구현 안됬습니다 ^__^\n";
        }
    };
};