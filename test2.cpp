#include "common.h"
#include <iostream>
#include <stdlib.h>
#include "custom_eval.h"
#include "lssvm.h"
#include <NTL/RR.h>
#include <random>
#include "hom_lssvm.h"





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

    EnDecoder Edcd(context);
    Encryptor enc(context);
    Decryptor dec(context);
    Bootstrapper btstr(eval);




    // Homomorphic matrix-vector multiplication
    cout << "Homomorphic matrix-vector multiplication!\n";
    long dim; 
    cout << "input dim: ";
    cin >> dim;
    long size; 
    cout << "input size: ";
    cin >> size;

    // plain matrix init
    cout << "plain matrix!\n";
    vector<vector<double>> matt(size,vector<double>(dim));
    for(int i=0; i<size; ++i) {
        for(int j=0; j<dim; ++j) {
            matt[i][j] = double(i*j-i);
            cout << matt[i][j] << " ";
        }
        cout << '\n';
    }
    cout << '\n';

    // plain vector init
    cout << "plain vector init!\n";
    vector<double> vec(size);
    for(int i=0; i<size; ++i) {
        cout << i << " ";
        vec[i] = (double)i;
    }
    cout << '\n';

    // Encrypt matrix
    cout << "Encrypt matrix!\n";
    vector<Ciphertext*> rows(size);
    for(int i=0; i<size; ++i) {
        rows[i] = new Ciphertext(context);
    }
    
    for(int i=0; i<size; ++i) {
        Message msg(logslot,0.0);
        for(int j=0; j<dim; ++j) {
            msg[j] = matt[j][(i+j)%dim];
            //printf("(%d,%d): ",j,(i+j)%size);
            //cout << msg[j].real() << " ";
        }
        //cout << '\n';
        enc.encrypt(msg,sk,*rows[i]);
    }

    HEmatrix A(context); 
    A.setHEmatrix(rows,size);

    // Encrypt vector
    cout << "Encrypt vector!\n";
    Ciphertext* v = new Ciphertext(context);
    Message msg(logslot,0.0);

    for(int i=0; i<size; ++i) {
        msg[i] = vec[i];
    }

    enc.encrypt(msg,sk,*v);

    HEvec b(context);
    b.setHEvec(v,size);

    cout << "init lssvm!\n";
    eval.init_lssvm(A,b,dim,size);

    
    return 0;
}