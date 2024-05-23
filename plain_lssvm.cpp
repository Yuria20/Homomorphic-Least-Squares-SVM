#include <iostream>
#include "lssvm.h"


int main() {

    /*
        X_ : vector<double>(size,vector<double(dim)) 
            반드시 데이터는 정규화 되어있어야 함.

        Y_ : vector<double>(size)
            Label은 1,-1 인 binary 형식이어야 함.

        dim : dimension of feature
        size : # of training sample
        constraint = 0.1 : Least Square SVM problem's constraint
        gamma = 0.1 : HyperParameter of rbf kernel.

        #############################
        gradient_descent 10만번 수행하고 1000번째마다 Loss := L2norm(b[n+1]-b[n]) 출력
        fit() 함수와 validation() 함수 참고
    */


    lssvm model(X_,Y_,dim,size,constraint,gamma);
    vector<double> pred = model.predict(x_,size);

    
    return 0;
}