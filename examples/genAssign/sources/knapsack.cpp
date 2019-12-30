#include "knapsack.h"
#include <cstring>
#include <cmath>
#include <new>
//#include <iostream>


namespace KNAPSACK {

double intKnapsack(int n, double *c, int* a, int b, int *x, double *dpMem)
{
	int beta, b0{0};
	double v,w, inf{-1.0e10}, zero{1.0e-10};
	double *F=(dpMem)? dpMem: new double[b+1];
	F[0]=w=0.0;
	for (beta=1; beta <= b; ++beta) {
		v=inf;
		for (int j=0; j < n; ++j) {
			if (beta >= a[j]) {
				if (v < F[beta-a[j]]+c[j])
					v=F[beta-a[j]]+c[j];
			}
		}
		if ((F[beta]=v) > w) {
			w=v;
			b0=beta;
		}
	}
// reverse step
	v=w; beta=b0;
	memset(x,0,n*sizeof(int));
	while (beta > 0) {
		for (int j{0}; j < n; ++j) {
			if (beta >= a[j])
				if (fabs(v - F[beta-a[j]] - c[j]) < zero) {
					x[j]++;
					v=F[beta-=a[j]];
					break;
				}
		}
	}
	if (!dpMem)
		delete[] F;
	return w;
} // end of intKnapsack

const int int_bit_size=sizeof(int)<<3;

inline int getBit(int k, int b, int *ipBt, int level_size)
{
	k=k*level_size + b / int_bit_size;
	return ipBt[k] & (1 << (b % int_bit_size));
}

inline void setBit(int k, int b, int *ipBt, int level_size)
{
	k=k*level_size + b / int_bit_size;
	ipBt[k]|=(1 << (b % int_bit_size));
} // end of setBit

void getMemForBinKnapsack(int n, int b, double* &dpF, int* &ipBt)
{
	int beta{b+1};
	int level_size{beta/int_bit_size};
	if (beta % int_bit_size)
		++level_size;
	dpF=new double[2*beta];
	ipBt=new int[n*level_size];
} // end of getMemForBinKnapsack

double binKnapsack(int n, double *c, int* a, int b, int *x, double *dpF, int *ipBt)
{
	double w, inf{-1.0e20};
	double *F1{dpF}, *F2;
	int beta{b+1};
	int level_size{beta/int_bit_size};
	F2=F1+beta;
	
	if (beta % int_bit_size)
		++level_size;
	memset(ipBt,0,n*level_size*sizeof(int)); // reset bit array
	F1[0]=0.0;
	for (beta=1; beta <= b; ++beta)
		F1[beta]=inf;
	F1[a[0]]=c[0];
	setBit(0,a[0],ipBt,level_size);
	for (int k, i{1}; i < n; ++i) {
		k=a[i]; w=c[i];
		for (beta=0; beta < k; beta++)
			F2[beta]=F1[beta];
		for (/*beta=k*/; beta <= b; ++beta) {
			if(F1[beta-k]+w > F1[beta]) {
				F2[beta]=F1[beta-k]+w;
				setBit(i,beta,ipBt,level_size);
			}
			else F2[beta]=F1[beta];
		}
		double *F{F1}; F1=F2; F2=F;
	}
// reverse step
	w=inf;
	for (int k{0}; k <= b; ++k) {
		if (w < F1[k])
			w=F1[beta=k];
	}
	for (int i{n-1}; i >= 0; --i) {
		if (getBit(i,beta,ipBt,level_size)) {
			x[i]=1;
			beta-=a[i];
		}
		else x[i]=0;
	}
	return w;
} // end of binKnapsack

};
