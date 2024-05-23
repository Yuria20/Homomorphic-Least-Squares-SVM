#pragma once

#include "common.h"
//#include"func.h"

using namespace std;
using namespace HEaaN;

class Poly {
public:
	double* coeff = 0;
	long deg, heap_l, heap_m, heaplen;
	double* chebcoeff = 0;
	Poly* poly_heap = 0;
	
	//Evaltree evaltree();
	Poly();
	Poly(long _deg);
	Poly(long _deg, double* _coeff, string tag);
	~Poly();

	void copy(Poly &poly);
	void set_polynomial(long _deg, double* _coeff, string tag);
	void showcoeff();
	void showchebcoeff();
	void power_to_cheb();
	void cheb_to_power();
	void set_zero_poly(Poly &ax);
	void set_zero_poly_special(Poly &ax);
	double evaluate(double value);
	Ciphertext evaluate(Ciphertext&ctxt);
	void construct_chebmat(double** &to_cheb);
	void gen_chebbasis(long deg,Poly& rtn);
	void construct_tree();
	double eval_heap(double &value);
	double traveling_tree(int idx,double &value,long &end);
	void add_poly(Poly &ax, Poly &bx, Poly &rtn);
	void sub_poly(Poly &ax, Poly &bx, Poly &rtn);
	void const_mul_poly(double c,Poly& rtn);
	void mul_poly(Poly &ax,Poly &bx,Poly &rtn);
	void power_poly(long deg, Poly& rtn);
	void divide_poly(Poly &px,Poly &tx,Poly &qx,Poly &rx);
};


