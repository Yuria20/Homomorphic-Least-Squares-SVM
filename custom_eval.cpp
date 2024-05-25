#include "custom_eval.h"
#include "poly.h"

void CustomEvaluator::HomPower(Ciphertext &ctxt, long deg, Ciphertext &rtn)
{   
    // Todo : binary multiplication optimization
    rtn = ctxt;
    for (int i = 1; i < deg; ++i)
    { 
        mult(rtn, ctxt, rtn);
    }
}

void CustomEvaluator::HomPolyEval(Ciphertext &ctxt, Poly &poly, Ciphertext &rtn)
{
    Ciphertext tmp(*zero);
    add(tmp, poly.coeff[0], tmp);

    for (int i = 1; i < poly.deg + 1; ++i)
    {
        if (poly.coeff[i] == 0.0)
            continue;

        Ciphertext val(*zero);

        // tmp = (ct)^n
        HomPower(ctxt, i, val);

        // tmp = a[i] * (ct)^n
        mult(val, poly.coeff[i], val);

        // rtn += tmp;
        add(val, tmp, tmp);
    }

    rtn = tmp;
}

void CustomEvaluator::set_exp(long exp_deg)
{   
    // init exp for some deg
    exp.deg = exp_deg;
    exp.set_zero_poly(exp);
    exp.coeff[0] = 1.0;
    for (int i = 1; i < exp_deg + 1; ++i)
    {
        exp.coeff[i] = exp.coeff[i - 1] / (double)i;
    }

    exp.power_to_cheb();
    exp.construct_tree();
}

void CustomEvaluator::trace(Ciphertext &ct,Ciphertext &rtn) {

    Ciphertext tmp(ct);
    Ciphertext tmp2(ct);

    //cout << logslot << '\n';
    for(int i=0; i<logslot; ++i) {
        //cout << (1<<i)<<"-leftRotate!\n";
        tmp2 = tmp;
        leftRotate(tmp2,(1<<i),tmp2);
        add(tmp,tmp2,tmp);
    }

    rtn = tmp;
}

void CustomEvaluator::trace(HEvec& vec, HEvec& rtn) {
    Ciphertext tmp = *vec.ct;
    rtn = vec;
    trace(tmp,tmp);
    *rtn.ct = tmp;
}

void CustomEvaluator::mat_vec_Multiplication(HEmatrix& mat, HEvec& vec, HEvec& rtn){
    if(mat.matsize != vec.dim) {
        cout << "size error!\n";
        return ;
    }

    if(mat.matsize < slot/2){
        Message hot(logslot,0.0);
        for(int i=0; i<mat.matsize; ++i) {
            hot[i] = 1.0;
        }

        Ciphertext tmp(*vec.ct); // (a b c d)
        Ciphertext tmp2(tmp);
        Ciphertext tmp3(*zero);
        mult(tmp,hot,tmp); // (a b 0 0 )
        rightRotate(tmp,mat.matsize,tmp2); // (0 0 a b)
        add(tmp,tmp2,tmp); // (a b a b)

        for(int i=0; i<mat.matsize; ++i) {
            mult(*mat.mat[i],tmp,tmp2);
            add(tmp3,tmp2,tmp3);
            leftRotate(tmp,1,tmp);
        }

        rtn = vec;
        *rtn.ct = tmp3;
    }
    else {

    }
};

void CustomEvaluator::gen_linear() {

    gen_kernel();

    Ciphertext y = *model->y_train->ct;

    long first = 0;
    long end = model->sample_size;

    for(int i=1; i<model->sample_size+1; ++i) {
        Ciphertext tmp(context);
        Message onehot(logslot,0.0);
        onehot[first] = 1.0;
        mult(y,onehot,tmp);
        leftRotate(tmp,first,tmp);
        add(*(model->linearMat->mat[i]),tmp,*(model->linearMat->mat[i]));

        Ciphertext tmp2(context);
        Message onehot2(logslot,0.0);
        onehot2[end-1] = 1.0;
        mult(y,onehot2,tmp2);
        rightRotate(tmp2,1,tmp2);
        add(*(model->linearMat->mat[i]),tmp2,*(model->linearMat->mat[i]));

        first++;
        end--;
    }
}

void CustomEvaluator::gen_kernel(){
    
    double per = 0.0;
    double num = (model->sample_size)*(model->sample_size);
    double delta = 1.0/num;
    auto start =std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    for(int i=0; i<model->sample_size; ++i) {
        for(int j=0; j<model->sample_size; ++j) {
            start =std::chrono::high_resolution_clock::now();

            Ciphertext tmp(context);
            sub(*(model->x_train->mat[i]),*(model->x_train->mat[j]),tmp);
            square(tmp,tmp);

            Message msg(logslot,0.0);
            for(int idx = 0; idx<model->feature_dim; ++idx) {
                msg[idx] = 1.0;
            }

            mult(tmp,msg,tmp);

            trace(tmp,tmp);

            negate(tmp,tmp);

            //multiplay Hyper Parameter gamma
            mult(tmp,model->gamma,tmp);

            Hom_exp(tmp,tmp);

            Message onehot(logslot,0.0);
            onehot[i] = 1.0;

            mult(tmp,onehot,tmp);
            
            mult(tmp,*(model->y_train->ct),tmp);

            Ciphertext Y(context);
            leftRotate(*(model->y_train->ct),j,Y);
            rightRotate(Y,i,Y);

            mult(tmp,Y,tmp);
            
            rightRotate(tmp,1,tmp);

            long idx = (i-j>0)?(model->sample_size+1-(i-j)):j-i;
            idx = idx%(model->sample_size+1);
            add(*(model->linearMat->mat[idx]),tmp,*(model->linearMat->mat[idx]));

            end = std::chrono::high_resolution_clock::now();
            num-=1.0;
            std::chrono::duration<double> duration = end - start;

            per+=delta;
            printf("kernel matrix... [%0.2f %%]\t Estimated remaining time:%f minutes\n",per*100, duration*60*num);
        }
    }

    Message I(logslot,1.0/model->constraint);
    I[0] = 0.0;
    add(*(model->linearMat->mat[0]),I,*(model->linearMat->mat[0]));

    cout << "kernel matrx done\n";
}

void CustomEvaluator::fit() {
    long epochs = 100;

    HEvec prev(context);
    long size = model->sample_size+1;

    auto start =std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    for(int i=0; i<epochs; ++i) {
        start = std::chrono::high_resolution_clock::now();

        Ciphertext tmp(*model->weight->ct);
        prev.setHEvec(&tmp, size);

        mat_vec_Multiplication(*model->linearMat,*model->weight,*model->weight);

        sub(*model->weight->ct,*model->tilda_one->ct,*model->weight->ct);

        mat_vec_Multiplication(*model->linearMat,*model->weight,*model->weight);

        mult(*model->weight->ct,model->learning_rate,*model->weight->ct);

        sub(*prev.ct,*model->weight->ct,*model->weight->ct);

        //cout << (*model->weight->ct).getLevel();    
        if((*model->weight->ct).getLevel()<7)
            btstr->bootstrapExtended(*model->weight->ct,*model->weight->ct);   


        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "(unit epochs)Execution time: " << duration.count() << " seconds" << std::endl;
    }
    
}
























void CustomEvaluator::gen_cheb_I(Ciphertext &ctxt, long heap_I)
{   
    // if there exist cheb_I, nothing
    if (cheb_I != 0)
        return;

    // memeory allocation
    long len = 1 << heap_I;
    cheb_I = new Ciphertext *[len];
    for (int i = 0; i < len; ++i)
    {
        cheb_I[i] = new Ciphertext(context);
    }

    // evaluate cheb[0] & cheb[1]
    *cheb_I[0] = *zero;
    add(*cheb_I[0], 1.0, *cheb_I[0]);
    *cheb_I[1] = Ciphertext(ctxt);

    // evaluate cheb[i] for i = 2,3, ... , 2 pow heap_I-1
    for (int i = 2; i < len; ++i)
    {
        mult(*cheb_I[i - 1], ctxt, *cheb_I[i]);
        mult(*cheb_I[i], 2.0, *cheb_I[i]);
        sub(*cheb_I[i], *cheb_I[i - 2], *cheb_I[i]);
    }
}

void CustomEvaluator::del_cheb_I(long heap_I){
    long len = 1 << heap_I;
    for(int i=0; i<len; ++i) {
        delete cheb_I[i];
        cheb_I[i] = nullptr;
    }
    delete [] cheb_I;
    cheb_I = nullptr;
}

void CustomEvaluator::gen_cheb_M(Ciphertext &ctxt, long heap_m, long heap_I)
{
    // if there not exist cheb_I, nothing
    if (cheb_I == 0)
        gen_cheb_I(ctxt,heap_I);

    // if there exist cheb_M, nothing
    if (cheb_M != 0)
        return;

    // memory allocation
    cheb_M = new Ciphertext *[heap_m - heap_I];
    for (int i = 0; i < heap_m - heap_I; ++i)
    {
        cheb_M[i] = new Ciphertext(context);
    }

    // evaluate cheb[I]
    long start = 1 << heap_I; 
    long end = 1 << heap_m;
    
    mult(*cheb_I[start - 1], ctxt, *cheb_M[0]);
    mult(*cheb_M[0], 2.0, *cheb_M[0]);
    sub(*cheb_M[0], *cheb_I[start - 2], *cheb_M[0]);
 
    // evaluate cheb[i] for i = 2^I, 2^I+1 , .... , 2^M-1
    for (int i = 1; i < heap_m - heap_I; ++i)
    {
        mult(*cheb_M[i - 1], *cheb_M[i - 1], *cheb_M[i]);
        mult(*cheb_M[i], 2.0, *cheb_M[i]);
        sub(*cheb_M[i], 1.0, *cheb_M[i]);
    }

}

void CustomEvaluator::del_cheb_M(long heap_I,long heap_m){
    long len = heap_m - heap_I;
    for(int i=0; i<len; ++i) {
        delete cheb_M[i];
        cheb_M[i] = nullptr;
    }
    delete [] cheb_M;
    cheb_M = nullptr;
}

void CustomEvaluator::BabystepEval(Ciphertext &ctxt, Poly &poly, Ciphertext &rtn)
{

    Ciphertext tmp(*zero);
    add(tmp, poly.chebcoeff[0], tmp);

    long len = poly.deg + 1;
    
    for (int i = 1; i < len; ++i){
        Ciphertext val(context);
        // var = a[i] * T[i]
        mult(*cheb_I[i], poly.chebcoeff[i], val);
        
        // rtn += tmp;
        add(val, tmp, tmp);
    }

    rtn = tmp;
}

void CustomEvaluator::eval_heap(Ciphertext &ctxt, Poly &poly, Ciphertext &rtn)
{
    if (poly.poly_heap == 0)
        poly.construct_tree();

    
    //ctheap.resize(poly.heaplen);
    //for(int i=0; i<ctheap.size(); ++i) {
        //ctheap[i] = new Ciphertext(context);
    //}
    //cout << "size: " << ctheap.size() << '\n';

    gen_cheb_I(ctxt, poly.heap_l);;
    gen_cheb_M(ctxt, poly.heap_m, poly.heap_l);

    //cout << "m:" << poly.heap_m << ",  l: "<< poly.heap_l << '\n';
    long end = (1 << (poly.heap_m - poly.heap_l + 1)) - 1;

    traveling_heap(ctxt, 1, rtn, end, poly);

    del_cheb_I(poly.heap_l);
    del_cheb_M(poly.heap_l,poly.heap_m);    
}

void CustomEvaluator::traveling_heap(Ciphertext &ctxt, long idx, Ciphertext &rtn, long end, Poly &poly)
{   //cout << idx << '\n';

    if (2 * idx > end)
    {   
        Ciphertext ttmp(rtn);
        BabystepEval(ctxt, poly.poly_heap[idx], ttmp); 
        //cout << "여기맞지? idx: "<< idx << '\n';
        //*ctheap[idx] = ttmp;
        rtn = ttmp;
        return;
    }

    Ciphertext tmp(rtn);
    Ciphertext tmp1(rtn), tmp2(rtn);
    traveling_heap(ctxt, 2 * idx, tmp1, end, poly);
    traveling_heap(ctxt, 2 * idx + 1, tmp2, end, poly);

    long index = floor(log2(2 * idx));
    mult(tmp1, *cheb_M[poly.heap_m - poly.heap_l - index], tmp);
    add(tmp, tmp2, tmp);

    rtn = tmp;
}

