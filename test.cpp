#include "poly.h"
#include <iostream>
#include <stdlib.h>
#include "custom_eval.h"
#include "lssvm.h"
#include <NTL/RR.h>
#include <random>

using namespace std;


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


    cout << "Homomorphic function Test Start!\n";
    double val = 2.0;
    Message z(logslot,val);
    Ciphertext ct(context);
    enc.encrypt(z,sk,ct);
    Ciphertext rtn(context);
    Message rtnmsg;
    int before;
    int after;
    long pownum;

    
    // Homomorphic power
    cout << "Homomorphic power evaluation test!\n";
    cout << "input degree: ";
    cin >> pownum;
    before = rtn.getLevel();
    eval.HomPower(ct,pownum,rtn);
    after = rtn.getLevel();
    dec.decrypt(rtn,sk,rtnmsg);
    printf("Homorphic Eval: %f,\t Real Eval: %f\n",rtnmsg[1].real(),pow(val,pownum));
    printf("%d-degree power evalution consume Depth %d!\n\n",pownum,before-after);

    // Homomprhic poly
    cout << "Homomorphic power evaluation test!\n";
    cout << "Test p(x) = 1 +2x +3x^2 + ... + nx^n-1\n";
    cout << "input degree: ";
    cin >> pownum;
    double *coeff = new double[pownum+1];
    for(int i=0; i<pownum+1; ++i) {
        coeff[i] = i+1.0;
    }
    Poly p(pownum,coeff,"power");
    rtn = ct;

    before = rtn.getLevel();
    eval.HomPolyEval(rtn,p,rtn);
    after = rtn.getLevel();

    dec.decrypt(rtn,sk,rtnmsg);
    printf("Homorphic Eval: %f,\t Real Eval: %f\n",rtnmsg[0].real(),p.evaluate(val));
    printf("%d-degree poly evalution consume Depth %d!\n\n",pownum,before-after);
    
    
    // Homomorphic exp
    cout << "Homomorphic exponential evaluation test!\n";
    cout << "input degree: ";
    cin >> pownum;
    eval.set_exp(pownum);

    rtn =  ct;

    before = rtn.getLevel();
    eval.Hom_exp(rtn,rtn);
    after = rtn.getLevel();

    dec.decrypt(rtn,sk,rtnmsg);
    printf("Homorphic Eval: %0.20f,\t Real Eval: %0.20f\n",rtnmsg[0].real(),exp(val));
    cout << "Error: "<< abs(rtnmsg[0].real()-exp(val)) << '\n';
    printf("%d-degree exp evalution consume Depth %d!\n\n",pownum,before-after);

    // Trace(Summation of all slots) operation
    cout << "Homomorphic trace evaluation test!\n";
    
    rtn = ct;
    
    before = rtn.getLevel();
    eval.trace(rtn,rtn);
    after = rtn.getLevel();

    dec.decrypt(rtn,sk,rtnmsg);
    printf("Homorphic trace: %f,\t Real Eval: %f\n",rtnmsg[0].real(),val*(1<<logslot));

    // Homomorphic matrix-vector multiplication
    cout << "Homomorphic matrix-vector multiplication!\n";
    long size; 
    cout << "input size: ";
    cin >> size;

    // plain matrix init
    cout << "plain matrix!\n";
    vector<vector<double>> matt(size,vector<double>(size));
    for(int i=0; i<size; ++i) {
        for(int j=0; j<size; ++j) {
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
        for(int j=0; j<size; ++j) {
            msg[j] = matt[j][(i+j)%size];
            printf("(%d,%d): ",j,(i+j)%size);
            cout << msg[j].real() << " ";
        }
        cout << '\n';
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

    eval.mat_vec_Multiplication(A,b,b);

    dec.decrypt(*b.ct,sk,msg);
    for(int i=0; i<size; ++i) {
        cout << msg[i].real() << " ";
    }
    cout << '\n';
    

    
    cout << "lssvm training!\n";
    long feature_dim,sample_size;
    cout << "feature_dim: "; cin >> feature_dim;
    cout << "sample_size: "; cin >> sample_size;
    cout << "random train_x & train_y gen!\n";

    vector<vector<double>> x_train(sample_size,vector<double>(feature_dim,0));
    vector<double> y_train(sample_size,0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.5,1.0);
    std::uniform_real_distribution<double> dis2(0.0,0.5);

    cout << "random x_train" << '\n';
    for(int i=0; i<sample_size; ++i) {
        for(int j=0; j<feature_dim; ++j) {
            if((i)%2) {
                x_train[i][j] = dis(gen);
                y_train[i] = 1.0;
            }
            else {
                x_train[i][j] = dis2(gen);
                y_train[i] = -1.0;
            }
            
            cout << x_train[i][j] << " ";
        }
        cout << '\n';
    }
    cout << '\n';

    cout << "random y_train" << '\n';
    //cout << y_train.size() << '\n';
    for(int j=0; j<sample_size; ++j) {
        //y_train[j] = (dis(gen)>0.5) ? 1.0 : -1.0;
        cout << y_train[j] << " ";
    }
    cout << '\n';
    
    lssvm model(x_train,y_train,feature_dim,sample_size,"rbf");

    cout << "test!!!\n";
    
    cout << "random x_test" << '\n';

    
    cout << "feature_dim: ";
    cin >> feature_dim;
    cout << "sample_size: ";
    cin >> sample_size;
    
    std::uniform_real_distribution<double> dis3(0.3,1.0);
    std::uniform_real_distribution<double> dis4(0.0,0.7);
    
    vector<vector<double>> x_test(sample_size,vector<double>(feature_dim,0));
    vector<double> y_test(sample_size,0);

    for(int i=0; i<sample_size; ++i) {
        for(int j=0; j<feature_dim; ++j) {
            if(i%2) {
                x_test[i][j] = dis3(gen);
                y_test[i] = 1.0;
            }
            else {
                x_test[i][j] = dis4(gen);
                y_test[i] = -1.0;
            }
            //cout << x_test[i][j] << " ";
        }
        //cout << '\n';
    }
    //cout << '\n';
    

    vector<double> pred = model.predict(x_test,sample_size);


    double acc = 0.0;
    for(int i=0; i<sample_size; ++i) {
        if(y_test[i] == pred[i]) {
            acc += 1.0;
            //cout << i <<"-result: " << y_test[i] << " " << pred[i] << '\n';
        }
    }

    //cout << sample_size;
    acc /= (double)(sample_size);
    
       
    cout << "accuracy: " << acc << '\n';
    //model.show_linear();

    return 0;
}