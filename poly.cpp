#include "poly.h"


Poly::Poly(){

}

Poly::Poly(long _deg) {
	deg = _deg;
	coeff = new double[deg + 1];
	chebcoeff = new double[deg + 1];
	for(int i = 0; i < deg + 1; i++) {
		coeff[i] = 0.0;
		chebcoeff[i] = 0.0;
	}
}

Poly::Poly(long _deg, double* _coeff, string tag) {
	deg = _deg;
	coeff = new double[deg + 1];
	chebcoeff = new double[deg + 1];
	if (tag == "power") {
		for(int i = 0; i < deg + 1; i++) {
			coeff[i] = _coeff[i];
		}
		power_to_cheb();
	}

	else if (tag == "cheb") {
		for(int i = 0; i < deg + 1; i++) {
			chebcoeff[i] = _coeff[i];
		}
		cheb_to_power();
	}
}

Poly::~Poly(){
	if(coeff) delete[] coeff;
	if(chebcoeff) delete[] chebcoeff;
	coeff = nullptr;
	chebcoeff = nullptr;
}

void Poly::copy(Poly &poly) {
	deg = poly.deg;
	coeff = new double[deg + 1];
	chebcoeff = new double[deg + 1];
	for(int i = 0; i < deg + 1; i++) {
		coeff[i] = poly.coeff[i];
		chebcoeff[i] = poly.chebcoeff[i];
	}

}
void Poly::set_polynomial(long _deg, double* _coeff, string tag) {
	deg = _deg;
	coeff = new double[deg + 1];
	chebcoeff = new double[deg + 1];
	if (tag == "power") {
		for(int i = 0; i < deg + 1; i++) {
			coeff[i] = _coeff[i];
		}
		power_to_cheb();
	}

	else if (tag == "cheb") {
		for(int i = 0; i < deg + 1; i++) {
			chebcoeff[i] = _coeff[i];
		}
		cheb_to_power();
	}
}

void Poly::showcoeff() {
	for(int i = 0; i < deg + 1; i++) {
		cout << "term " << i << " : " << coeff[i] << endl;
	}
	cout << endl;
}

void Poly::showchebcoeff() {
	for(int i = 0; i < deg + 1; i++) {
		cout << "chebterm " << i << " : " << chebcoeff[i] << endl;
	}
	cout << endl;
}

void Poly::power_to_cheb() {
	// construct cheb_basis
	double** to_cheb=nullptr;
	construct_chebmat(to_cheb);

	// copy coeff
	double* tmp = new double[deg+1];
	for(int i=0; i<deg+1; ++i) {
		tmp[i] = coeff[i];
	}

	// gaussian elimination 
	for(int i=deg; i>=0; --i) {
		chebcoeff[i] = tmp[i]/to_cheb[i][i];
		for(int k=0; k<i; ++k) {
			tmp[k]-=(to_cheb[i][k]*chebcoeff[i]);
		}
	}

	// delete dynamic memory
	for(int i=0; i<deg+1; ++i) {
		delete [] to_cheb[i];
		to_cheb[i] = nullptr;
	}
	delete [] to_cheb;
	to_cheb = nullptr;
	delete [] tmp;
	tmp = nullptr;
}

void Poly::cheb_to_power() {
	// construct cheb_basis
	double** to_cheb=nullptr; 
	construct_chebmat(to_cheb);

	// matrix-vector multiplication
	for(int i=0; i<deg+1; ++i) {
		double tmp = 0.0;
		for(int j=0; j<deg+1; ++j) {
			tmp+=(coeff[j]*to_cheb[j][i]);
		}
		chebcoeff[i] = tmp;
	}

	// delete
	for(int i=0; i<deg+1; ++i) {
		delete [] to_cheb[i];
		to_cheb[i]=nullptr;
	}
	delete [] to_cheb;
	to_cheb = nullptr;
	
}

void Poly::set_zero_poly(Poly &ax) {

	ax.coeff = new double[ax.deg+1];
	ax.chebcoeff = new double[ax.deg+1];

	ax.coeff[0]=0.0;

	for(int i=0; i<ax.deg+1; ++i) {
		ax.coeff[i]=0.0;
	}

}

void Poly::set_zero_poly_special(Poly &ax) {
	ax.deg = 0;

	ax.coeff = new double[ax.deg+1];
	ax.chebcoeff = new double[ax.deg+1];

	ax.coeff[0]=0.0;

	for(int i=0; i<ax.deg+1; ++i) {
		ax.coeff[i]=0.0;
	}

}

double Poly::evaluate(double value) {
	double rtn = 0.0, term = 1.0;
	for(int i = 0; i <= deg; i++) {
		rtn += coeff[i] * term;
		term *= value;
	}
	
	return rtn;
}

Ciphertext Poly::evaluate(Ciphertext &ctxt){
	for(int i=0; i<deg; ++i) {
	}
}

void Poly::construct_chebmat(double** &to_cheb) {
		
	// memory allocation
	to_cheb = new double*[deg+1];
	for(int i=0; i<deg+1; ++i) {
		to_cheb[i] = new double[deg+1];
	}

	for(int i=0; i<deg+1; ++i) {
		for(int j=0; j<deg+1; ++j) {
			to_cheb[i][j]=0;
		}
	}

	// init T[0]=1, T[1]=x
	to_cheb[0][0] = 1.0;
	if(deg==0)
		return;

	to_cheb[1][1] = 1.0;
	if(deg==1)
		return;

	// recursively compute T[n+2] = 2xT[n+1]-T[n]
	for(int i=2; i<deg+1; ++i) {
		for(int k = 1; k<deg+1; ++k) {
			to_cheb[i][k] += 2*to_cheb[i-1][k-1];
		}
		for(int k = 0; k<deg+1; ++k) {
			to_cheb[i][k] -= to_cheb[i-2][k] ;
		}
	}
}

void Poly::gen_chebbasis(long deg,Poly& rtn) {
	Poly tmp(deg);
	double** to_cheb=nullptr;
	tmp.construct_chebmat(to_cheb);

	for(int i=0; i<tmp.deg+1; ++i) {
		tmp.coeff[i] = to_cheb[deg][i];
	}

	rtn.copy(tmp);

	// delete
	for(int i=0; i<deg+1; ++i) {
		delete [] to_cheb[i];
		to_cheb[i]=nullptr;
	}
	delete [] to_cheb;
	to_cheb=nullptr;
}

void Poly::construct_tree(){
	// init heap 
	heap_m= ceil(log2(deg+1));
	heap_l = floor((double)heap_m/2);
	heaplen = (1<<heap_m)-1;

	poly_heap = new Poly[heaplen+1];
	poly_heap[1].copy(*this);
	//poly_heap[1].showcoeff();

	// construct cheb_mat
	double** to_cheb = nullptr; 
	construct_chebmat(to_cheb);

	// construct evaluation_tree
	long chebdeg = 1 << (heap_m-1);
	//cout << "chebdeg: " << chebdeg << '\n';
	//cout << "heaplen: " << heaplen << '\n';

	for(int i=1; i<((heaplen)>>1)+1; ++i) {
		
		Poly cheb(chebdeg,to_cheb[chebdeg],"power");
		//cout << "i:" << i << '\n';

		if(floor(log2(i))!=floor(log2(i+1))) {
			chebdeg = chebdeg >> 1;
		}

	
		if(poly_heap[i].deg >= cheb.deg) {
			divide_poly(poly_heap[i],cheb,poly_heap[2*i],poly_heap[2*i+1]);
		}
		else {
			//double cof[] = {0.0};
			//poly_heap[2*i] = Poly(0,cof,"power");
			poly_heap[2*i].set_zero_poly_special(poly_heap[2*i]);
			poly_heap[2*i+1].copy(poly_heap[i]);
		}

		//poly_heap[i].power_to_cheb();
	}

    for(int i=1; i<heaplen+1; ++i) {
        poly_heap[i].power_to_cheb();
    }
	
	// delete cheb_mat
	for(int i=0; i<deg+1; ++i) {
		delete [] to_cheb[i];
		to_cheb[i]=nullptr;
	}
	delete [] to_cheb;
	to_cheb=nullptr;

}

void Poly::add_poly(Poly &ax, Poly &bx, Poly &rtn) {

	long degree = ax.deg > bx.deg ? ax.deg : bx.deg;
	rtn.deg = degree;
	Poly tmp(degree);

	for(int i=0; i<ax.deg+1; ++i) {
		tmp.coeff[i]+=ax.coeff[i];
	}

	for(int i=0; i<bx.deg+1; ++i) {
		tmp.coeff[i]+=bx.coeff[i];
	}

	rtn.copy(tmp);
}

void Poly::sub_poly(Poly &ax, Poly &bx, Poly &rtn) {
	long degree = ax.deg > bx.deg ? ax.deg : bx.deg;
	Poly tmp(degree);

	for(int i=0; i<ax.deg+1; ++i) {
		tmp.coeff[i]+=ax.coeff[i];
	}

	for(int i=0; i<bx.deg+1; ++i) {
		tmp.coeff[i]-=bx.coeff[i];
	}

	long leading = degree;
	while(tmp.coeff[leading]==0 && leading!=0) {
		leading--;
	}

	tmp.deg = leading;
	rtn.copy(tmp);
}

void Poly::const_mul_poly(double c,Poly& rtn) {
	for(int i=0; i<rtn.deg+1;++i) {
		rtn.coeff[i]*=c;
	}
}

void Poly::mul_poly(Poly &ax,Poly &bx,Poly &rtn) {
	long degree = ax.deg+bx.deg;
	Poly tmp(degree);

	for(int i=0; i<ax.deg+1; ++i) {
		for(int j=0; j<bx.deg+1; ++j) {
			tmp.coeff[j+i] +=ax.coeff[i]*bx.coeff[j];
		}
	}

	rtn.copy(tmp);
}

void Poly::power_poly(long deg, Poly& rtn) {
	Poly power(1);
	power.coeff[1] = 1.0;
		
	Poly tmp(0);
	tmp.coeff[0]=1.0;

	for(int i=0; i<deg; ++i) {
		power.mul_poly(tmp,power,tmp);
	}

	rtn.copy(tmp);
}

void Poly::divide_poly(Poly &px,Poly &tx,Poly &qx,Poly &rx) {

	if(px.deg < tx.deg) {
		cout << "Error exeption : degree error!";
		return;
	}

	long qx_deg = px.deg-tx.deg;
	Poly q_tmp(qx_deg);
	q_tmp.set_zero_poly(q_tmp);

	long rx_deg = px.deg;
	Poly r_tmp; 
	r_tmp.copy(px);
	
	while(r_tmp.deg>=tx.deg) { 
	
		Poly tx_tmp;
		tx_tmp.copy(tx);
		
		double leading = r_tmp.coeff[r_tmp.deg]/tx_tmp.coeff[tx_tmp.deg];
		q_tmp.coeff[qx_deg] = leading;
		qx_deg--;	

		tx_tmp.const_mul_poly(leading,tx_tmp);

		Poly power;
		power.power_poly(r_tmp.deg-tx_tmp.deg, power);

		tx_tmp.mul_poly(tx_tmp,power,tx_tmp);
		r_tmp.sub_poly(r_tmp,tx_tmp,r_tmp);
	}

	qx.copy(q_tmp);
	rx.copy(r_tmp);
}

double Poly::eval_heap(double& value) {

	if (poly_heap == 0)
        construct_tree();

    //cout << "m:" << poly.heap_m << ",  l: "<< poly.heap_l << '\n';
    long end = (1 << (heap_m - heap_l + 1)) - 1;

    double ret = traveling_tree(1,value,end);

}

double Poly::traveling_tree(int idx, double&value, long& end) {
    //cout << "idx: " << idx << '\n'; 
    //poly_heap[idx].showchebcoeff();

	if(2*idx > end) {
		//poly_heap[idx].showchebcoeff();
		cout<< "poly_heap["<<idx<<"]: " << poly_heap[idx].evaluate(value) << '\n';
		return poly_heap[idx].evaluate(value);
	}

	long index = 1 << (int)floor(log2(heaplen)-log2(idx)) ;

	Poly tx(index);
	tx.gen_chebbasis(tx.deg,tx);
	//tx.showcoeff();

	double t1 = traveling_tree(2*idx,value, end);
	double t2 = traveling_tree(2*idx+1,value,end);
    //printf("%f*%f+%f\n",t1,tx.evaluate(value),t2);

	return t1*tx.evaluate(value)+t2;
}