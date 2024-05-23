#pragma once

#include "common.h"
#include "poly.h"
#include "HEmatrix.h"
#include "hom_lssvm.h"

class CustomEvaluator : public HomEvaluator {
public:
    Context context;
    long logslot;
    long slot;
    Ciphertext* zero;
    hom_lssvm model;

    /* Only Debugging */
    //vector<Ciphertext*> ctheap;

    // Constructor
    CustomEvaluator(const Context &context, const std::string &key_dir_path):HomEvaluator(context,key_dir_path){};
    CustomEvaluator(const Context &context, const KeyPack &pack):HomEvaluator(context,pack) {
        this->context = context;
        this->logslot = getLogFullSlots(context);
        this->slot = 1 << logslot;
        this->zero = new Ciphertext(this->context);
        sub(*zero,*zero,*zero);
    };

    // Polynomial Evalutation
    void HomPower(Ciphertext &ctxt, long deg, Ciphertext &rtn);
    void HomPolyEval(Ciphertext &ctxt, Poly &poly, Ciphertext &rtn);
    void HomPolyChebEval(Ciphertext &ctxt, Poly &poly, Ciphertext &rtn);

    void set_exp(long exp_deg);
    void Hom_exp(Ciphertext &ctxt, Ciphertext &rtn){
        eval_heap(ctxt, exp, rtn);
    }

    // Matrix-Vector Evaluation
    void trace(Ciphertext &ct,Ciphertext &rtn);
    void trace(HEvec& vec, HEvec& rtn);
    void mat_vec_Multiplication(HEmatrix& mat, HEvec& vec, HEvec& rtn);

    // lssvm operation
    void init_lssvm(HEmatrix& train_x, HEvec& train_y,long feature_dim, long sample_size) {
        cout << "fuck\n";
        
        (*model.x_train).mat.resize(sample_size);

        cout << "fuck\n";
        
        for(int i=0; i<sample_size; ++i) {
            (*model.x_train).mat[i] = new Ciphertext(context);
        }

        (*model.y_train).ct = new Ciphertext(context);


        *model.x_train = train_x;
        *model.y_train = train_y;

        model.feature_dim = feature_dim;
        model.sample_size = sample_size;

        
        // initial b = (0.1)
        Ciphertext b(*zero);
        Message msg(logslot,0.1);
        mult(b,msg,b);
        *((*model.weight).ct) = b;

        // tilda one init
        Ciphertext one(*zero);
        Message msg2(logslot,1.0);
        msg2[0] = 0.0;
        mult(one,msg2,one);
        *((*model.tilda_one).ct) = one;

        ((*model.linearMat).mat).resize(sample_size+1);
        for(int i=0; i<sample_size+1; ++i) {
            (*model.linearMat).mat[i] = new Ciphertext(context);
        }
    }

    void gen_linear();
    void gen_kernel();

private:
    Poly exp;
    
    Ciphertext **cheb_I = 0;
    Ciphertext **cheb_M = 0;

    void gen_cheb_I(Ciphertext &ctxt, long heap_I);
    void gen_cheb_M(Ciphertext &ctxt, long heap_m, long heap_I);
    void del_cheb_I(long heap_I);
    void del_cheb_M(long heap_I,long heap_m);
    void BabystepEval(Ciphertext &ctxt, Poly &poly, Ciphertext &rtn);
    void eval_heap(Ciphertext &ctxt, Poly &poly, Ciphertext &rtn);
    void traveling_heap(Ciphertext &ctxt, long idx, Ciphertext &rtn, long end, Poly &poly);
};