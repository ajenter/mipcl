#include <limits>
#include <fstream>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>
#include <Sort.h>
#include <except.h>
#include "cutstock.h"

CCutStock::CCutStock(const char* name): m_dpF(0), CMIP(name)
{
}

#ifndef __ONE_THREAD_
CCutStock::CCutStock(const CCutStock &other, int thread): CMIP(other,thread)
{
}

CMIP* CCutStock::clone(const CMIP *pMip, int thread)
{
	return static_cast<CMIP*>(new CCutStock(*static_cast<CCutStock*>(const_cast<CMIP*>(pMip)),thread));
}
#endif

CCutStock::~CCutStock()
{
#ifndef __ONE_THREAD_
	if (!m_iThread)
#endif
	if (m_dpF)
		delete[] m_dpF;
}

void CCutStock::init(int m, int *l, int *q, int L)
{
	int sz, *b, *S;
	m_iFinalTypeNum=m;
	m_ipFinalLength=l; 
	m_ipFinalNum=q;
	m_iRawLength=L;
	if (!(m_dpF = new double[L+1])) {
		throw new CMemoryException("CCutStock::init()");
	}
	openMatrix(m,m,m*m,false,true,m,3*m,3*m*m);
	S=reinterpret_cast<int*>(m_dpW);
	for (int i=0; i < m; ++i) {
		addCtr(S[i]=i,0,static_cast<double>(q[i]),CLP::INF);
	}

	memcpy(b=S+m,q,m*sizeof(int));
	SORT::decSortInt(m,S,l);
	for (int j{0}, i{m}; i; ) {
		newColumn(i,S,b,sz);
		addColumn(j++,VAR_INT,-1.0,0.0,VAR_INF,sz,m_dpArray,m_ipArray);
	}
	
	preprocOff();
	setObjSense(false);
//	setScaling(SCL_NO);
	closeMatrix();
} // end of CCutStock::init()

void CCutStock::newColumn(int &k, int *S, int *b, int &sz)
{
	int q, t, t0{std::numeric_limits<int>::max()},
		m{m_iFinalTypeNum}, r{m_iRawLength};
	double *dpVal{m_dpArray};
	int	*w{m_ipFinalLength}, *ipRow{m_ipArray};
	int *ipNum{reinterpret_cast<int*>(dpVal+m)};

	memset(ipNum,0,m*sizeof(int));
	sz=0;
	for (int i{0}; i < k; ++i) {
		int j{S[i]};
		if (r >= w[j]) {
			t=b[j]/(q=r/w[j]);
			if (b[j] % q)
				++t;
			if (t0 > t)
				t0=t;
			dpVal[sz]=static_cast<double>(ipNum[j]=q);
			ipRow[sz++]=j;
			if (!(r%=w[j]))
				break;
		}
	}

	q=0;
	for (int i{0}; i < k; ++i) {
		int j{S[i]};
		if ((b[j]-=t0*ipNum[j]) > 0)
			S[q++]=j;
	}
	k=q;
} // end of CCutStock::newColumn()

bool CCutStock::generateColumns(int m, const tagHANDLE* ipRowHd, const double* dpY)
{
	if (getCurrentNodeHeight())
		return false;
	double *c{m_dpArray+m};
	int *x{m_ipArray+m};
	bool flag{false};
	for (int i{0}; i < m; ++i)
		c[i]=-dpY[i];
	if (intKnapsack(m,c,m_ipFinalLength,m_iRawLength,x,m_dpF) > 1.0+getVarTol()) {
		int sz{0};
		for (int i{0}; i < m; ++i) {
			if (x[i] > 0) {
				m_dpArray[sz]=(double)x[i];
				m_ipArray[sz++]=i;
			}
		}
		addNewColumn(-1,VAR_INT,1.0,0.0,VAR_INF,
                     sz,m_dpArray,m_ipArray,
                     false,false,0,true);
		flag=true;
	}
	return flag;
} // end of CCutStock::generateColumns()

double CCutStock::intKnapsack(int n, double *c, int* a, int b, int *x, double *dpMem)
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
					++x[j];
					v=F[beta-=a[j]];
					break;
				}
		}
	}
	if (!dpMem)
		delete[] F;
	return w;
} // end of CCutStock::intKnapsack()

void CCutStock::printSolution(const char* name)
{
	char str[128];
	if (name)
		strcpy(str,name);
	else {
		getProblemName(str);
		strcat(str,".sol");
	}
	std::ofstream fout(str);
	if (!fout.is_open()) {
		throw new CFileException("CCutStock::printSolution",str);
	}
	fout << static_cast<int>(floor(getObjVal()+0.5)) << " rolls\n";
	fout << "Patterns:\n";
	double zero{std::numeric_limits<double>::epsilon()};
	double *x, *dpVal{m_dpArray};
	int *ipHd, *ipRow{m_ipArray};
	int n,sz;

	n=getSolution(x,ipHd);
	for (int j{0}; j < n; ++j) {
		if (x[j] > zero) {
			fout << x[j] << ":";
			sz=getColumn(j,dpVal,ipRow,false);
			for (int i=0; i < sz; ++i) {
				fout << " " << static_cast<int>(floor(dpVal[i]+0.5)) << " of " << ipRow[i]+1 << ";";
			}
			fout << std::endl;
		}
	}
	fout.close();
} // end of CCutStock::printSolution()

