#pragma once
#include "common.h"
#include "HEmatrix.h"
#include "custom_eval.h"

class hom_lssvm {
public:
    Context context;
    long feature_dim;
    long sample_size;
    HEmatrix* x_train;
    HEvec* y_train;
    double learning_rate = 0.0001;
    double constraint = 1.0;
    double gamma = 0.1;
    HEvec* weight;
    HEvec* tilda_one;
    HEmatrix* linearMat; // (sample_size+1) x (sample_size+1)
    
    hom_lssvm(Context context) {
        this->context = context;
        x_train = new HEmatrix(context);
        y_train = new HEvec(context);
        weight = new HEvec(context);
        tilda_one = new HEvec(context);
        linearMat = new HEmatrix(context);
    };
};