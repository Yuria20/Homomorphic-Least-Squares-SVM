#include "common.h"
#include <iostream>
#include <stdlib.h>
#include "custom_eval.h"
#include "lssvm.h"
#include <NTL/RR.h>
#include <random>
#include "hom_lssvm.h"

void GenRandomData(vector<vector<double>> &x_train, vector<double> &y_train,long size,long dim) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.5,1.0);
    std::uniform_real_distribution<double> dis2(0.0,0.5);

    bool isprint = false;

    cout << "gen random training data..." << '\n';
    for(int i=0; i<size; ++i) {
        for(int j=0; j<dim; ++j) {
            if((i)%2) {
                x_train[i][j] = dis(gen);
                y_train[i] = 1.0;
            }
            else {
                x_train[i][j] = dis2(gen);
                y_train[i] = -1.0;
            }
        }
    }

    if(isprint) {
        cout << "X_train:\n";
        for(int i=0; i<size; ++i) {
            for(int j=0; j<dim; ++j) {
                cout << x_train[i][j] << " ";
            }
            cout << '\n';
        }
        cout << "Y_train:\n";
        for(int i=0; i<size; ++i){
            cout << y_train[i] << " ";
        }
        cout << '\n';
    }

   
}

int main() {

    // Generate pre-determined Context where FVa
    HEaaN::ParameterPreset params = HEaaN::ParameterPreset::FVa;
    HEaaN::Context context = makeContext(params);
    printf("parameter setting complete!\n");

    long logslot = HEaaN::getLogFullSlots(context);
    long slot = 1 << logslot;
    printf("logslot: %ld,\tnumber of slot: %ld\n",logslot,slot);

    // Secret Key Gen
    HEaaN::SecretKey sk = HEaaN::SecretKey(context);
    printf("secret key generating complete!\n");

    // Public, Multiplication, Conjugation, Rotation Key Gen
    HEaaN::KeyGenerator keygen = HEaaN::KeyGenerator(context,sk);
    keygen.genEncryptionKey();
    printf("encryption key generating complete!\n");
    keygen.genMultiplicationKey();
    printf("mutliplication key generating complete!\n");
    keygen.genConjugationKey();
    printf("conjugation key generating complete!\n");
    keygen.genRotationKeyBundle();
    printf("rotation key generating complete!\n\n");
    keygen.genRotKeysForBootstrap(16);

    // Init Evaluator, Encoder, Decoder, Encryptor, Decryptor
    //HomEvaluator eval(context,keygen.getKeyPack());
    KeyPack keys = keygen.getKeyPack();
    CustomEvaluator eval(context,keygen.getKeyPack());
    eval.set_exp(15);

    EnDecoder Edcd(context);
    Encryptor enc(context);
    Decryptor dec(context);
    Bootstrapper btstr(eval);
    eval.btstr = &btstr;
    

    /* Homomorphic matrix-vector multiplication */ 

    cout << "Homomorphic Least Square SVM test!\n";
    long dim,size; 
    dim = size = 1000;
    //cout << "input dim: "; //cin >> dim;
    //cout << "input size: "; //cin >> size;

    vector<vector<double>> x_train(size,vector<double>(dim,0));
    vector<double> y_train(size,0);

    GenRandomData(x_train,y_train,size,dim);

    cout << "Encrypt matrix...\n";
    vector<Ciphertext*> X_(size);
    for(int i=0; i<size; ++i) {
        X_[i] = new Ciphertext(context);
    }
    
    for(int i=0; i<size; ++i) {
        Message msg(logslot,0.0);
        for(int j=0; j<dim; ++j) {
            msg[j] = x_train[i][j];
        }
        enc.encrypt(msg,sk,*X_[i]);
    }

    cout << "Encrypt vector...\n";
    Ciphertext* Y_ = new Ciphertext(context);
    Message msg(logslot,0.0);
    for(int i=0; i<size; ++i) {
        msg[i] = y_train[i];
    }
    enc.encrypt(msg,sk,*Y_);
    
    /*debug*/
    Ciphertext cc(context);
    Message mm(logslot);

    lssvm plain_model(x_train,y_train,dim,size);
    cout << "PlainText\n";
    //plain_model.show_linear();

    cout << "lssvm model init\n";
    hom_lssvm model(context);
    eval.set_lssvm(&model);
    eval.init_lssvm(X_,Y_,dim,size);

    cout << "gen_linear...\n";

    auto start = std::chrono::high_resolution_clock::now();
    eval.gen_linear();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::ratio<60>> duration = end - start;
    std::cout << "(kernel setting)Execution time: " << duration.count() << " minutes" << std::endl;

    eval.fit();

    //debug 용도 코드
    cout << "PlainText\n";
    plain_model.show_weight();
    cout << "CipherText\n";
    dec.decrypt(*(eval.model->weight->ct),sk,mm);

    for(int i=0; i<size+1; ++i) {
        cout << mm[i].real() << " ";
    }

    return 0;
}