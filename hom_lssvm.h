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
    double learning_rate = 0.00000001;
    double constraint = 0.1;
    HEvec* weight;
    HEvec* tilda_one;
    HEmatrix* linearMat; // (sample_size+1) x (sample_size+1)
    //string tag = "linear";
};