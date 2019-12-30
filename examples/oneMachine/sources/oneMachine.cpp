#include <limits>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstring>
#include <except.h>
#include <Sort.h>
#include "oneMachine.h"

COneMachine::COneMachine(const char* name, int n, double *c, int *p, int *r, int* d): CMIP(name)
{
	int i, sz,dmax=-1, rmin=std::numeric_limits<int>::max();
	m_iJobNum=n;
	m_ipP=p; m_ipR=r; m_ipD=d; m_dpCost=c;
	openMatrix(6*n,2*n,n+(n*n)/4);
    for (i=sz=0; i < n; i++) {
        addVar(i,VAR_BIN,c[i],0.0,1.0);
        if (p[i] > 0) {
            m_dpArray[sz]=p[i];
            m_ipArray[sz++]=i;
        }
        if (r[i] < rmin)
            rmin=r[i];
        if (d[i] > dmax)
            dmax=d[i];
    }
    m_iDmax=dmax;
    for (i=0; i < n; i++) {
    	setVarPriority(addVar(i+n,VAR_INT,0.0,r[i],d[i]-p[i]),VAR_PRI_MIN);
    }
    addRow(0,CTR_KNAPSACK | CTR_INT,-INF,dmax-rmin,sz,m_dpArray,m_ipArray);

    preprocOff();
    setScaling(SCL_NO);
    closeMatrix();
} // end of COneMachine::COneMachine()

#ifndef __ONE_THREAD_
COneMachine::COneMachine(const COneMachine &other, int thread): CMIP(other,thread)
{
	m_iJobNum=other.m_iJobNum;
	m_ipP=other.m_ipP;
	m_ipR=other.m_ipR;
	m_ipD=other.m_ipD;
	m_dpCost=other.m_dpCost;
	m_iDmax=other.m_iDmax;
}

CMIP* COneMachine::clone(const CMIP *pMip, int thread)
{
	return static_cast<CMIP*>(new COneMachine(*static_cast<COneMachine*>(const_cast<CMIP*>(pMip)),thread));
}
#endif

bool COneMachine::isFeasible(int varNum, const double* dpX, const tagHANDLE* ipColHd)
{
    int i,n,j1,j2, *ipStart, *ipTsk;
    n=m_iJobNum;
    ipTsk=(ipStart=m_ipArray)+n;
    for (i=j1=0; i < n; i++) {
        if (dpX[i] > 0.5)
            ipStart[ipTsk[j1++]=i]=(int)floor(dpX[n+i]+0.5);
    }
    SORT::incSortInt(n=j1,ipTsk,ipStart);
    j1=ipTsk[0];
    for (i=1; i < n; i++) {
        j2=ipTsk[i];
        if (ipStart[j1]+m_ipP[j1] > ipStart[j2])
            return false;
        j1=j2;
    }
    return true;
} // end of COneMachine::isFeasible()

int COneMachine::startBranching(int nodeHeight)
{
    int k;
    m_iJob1=-1;
    if (!(k=CMIP::startBranching(nodeHeight))) {
        double *dpX;
        int *ipHd, *ipStart, *ipTsk;
        int n,j1,j2;
        n=m_iJobNum;
        ipTsk=(ipStart=m_ipArray)+n;
        CLP::getSolution(dpX,ipHd);
        for (int i=k=0; i < n; ++i) {
            if (dpX[i] > 0.5)
                ipStart[ipTsk[k++]=i]=(int)floor(dpX[n+i]+0.5);
        }
        SORT::incSortInt(k,ipTsk,ipStart);
        j1=ipTsk[0];
        n=k; k=0;
        for (int i=1; i < n; ++i) {
            j2=ipTsk[i];
            if (ipStart[j1]+m_ipP[j1] > ipStart[j2]) {
                k=3;
                m_iJob1=j1;
                m_iJob2=j2;
                break;
            }
            j1=j2;
        }
    }
    return k;
} // end of COneMachine::startBranching()

bool COneMachine::updateBranch(int k)
{
    bool flag;
    if (m_iJob1 < 0) {
        flag=CMIP::updateBranch(k);
    }
    else {
    	double dpVal[2];
        int hd, j1, j2, ipCol[2];
        if (k) {
            j1=m_iJob1;
            j2=m_iJob2;
        }
        else {
            j1=m_iJob2;
            j2=m_iJob1;
        }
        hd=(j1 | (j2 << (4*sizeof(int)))) << 1;
        if (k < 2) {
            dpVal[0]=1.0; dpVal[1]=-1.0;
            ipCol[0]=m_iJobNum+j2;
            ipCol[1]=m_iJobNum+j1;
            addCut(hd,CTR_INT|CTR_LOCAL|CTR_ATTACHED,m_ipP[j1],INF,2,dpVal,ipCol);
            setVarLoBound(j1,1.0);
            setVarLoBound(j2,1.0);
        }
        else {
            dpVal[0]=dpVal[1]=1.0;
            ipCol[0]=j1; ipCol[1]=j2;
            addCut(hd+1,CTR_INT|CTR_PACKING|CTR_LOCAL|CTR_ATTACHED,-INF,1.0,2,dpVal,ipCol);
        }
        flag=true;
    }
    return flag;
} // end of COneMachine::updateBranch()

bool COneMachine::getRow(tagHANDLE hd,
        int varNum, const tagHANDLE* ipColHd,
        int& type, double& b1, double& b2,
        int& sz, double* dpVal, int* ipCol, bool &scaled)
{
    int k,j1, j2;
    k=hd & 0x1;
    hd>>=1;
    j2=hd >> (4*sizeof(int));
    j1=hd ^ (j2 << (4*sizeof(int)));

    type=CTR_INT | CTR_LOCAL | CTR_ATTACHED;
    sz=2;
    if (k == 1) {
        type|=CTR_PACKING;
        b1=-INF;
        b2=1.0;
        dpVal[0]=dpVal[1]=1.0;
        ipCol[0]=j1;
        ipCol[1]=j2;
    }
    else {
        b2=INF;
        b1=m_ipP[j1];
        dpVal[0]=1.0;
        dpVal[1]=-1.0;
        ipCol[0]=m_iJobNum+j2;
        ipCol[1]=m_iJobNum+j1;
    }
    scaled=true;
    return true;
} // end of COneMachine::getRow()

bool COneMachine::separate(int varNum,
        const double* dpX, const tagHANDLE* colHd, bool genFlag)
{
	double w;
	int *ipR, *ipD, *ipT1, *ipT2;
	int k1,k2,d,delta,s1,s2,sz,type,t1,t2,cutNum=0,n=m_iJobNum;

	ipT1=reinterpret_cast<int*>(m_dpFd); ipT2=ipT1+n;
	ipR=reinterpret_cast<int*>(m_dpFb); ipD=ipR+n;

	for (int i=sz=0; i < n; ++i) {
		ipR[i]=static_cast<int>(floor(getVarLoBound(i+n)+0.5));
		ipD[i]=static_cast<int>(floor(getVarUpBound(i+n)+0.5))+m_ipP[i];
		if (dpX[i] > 0.1) {
			ipT1[sz]=ipR[i];
			ipT2[sz++]=-ipD[i];
		}
	}
	std::sort(ipT1,ipT1+sz);
	std::sort(ipT2,ipT2+sz);
	t1=ipT1[0]; t2=(ipT2[0]=-ipT2[0]);
	for (int i=k1=k2=1; i < sz; ++i) {
		if (t1 < ipT1[i])
			t1=ipT1[k1++]=ipT1[i];
		ipT2[i]=-ipT2[i];
		if (t2 > ipT2[i])
			t2=ipT2[k2++]=ipT2[i];
	}
	for (int i=0; i < k1; ++i) {
		t1=ipT1[i];
		for (int j=0; j < k2; ++j) {
			t2=ipT2[j];
			if ((delta=t2-t1) <= 0)
				break;
			w=0.0;
			type=CTR_INT | CTR_KNAPSACK;
			for (int k=sz=0; k < n; ++k) {
				s1=s2=m_ipP[k];
				if (ipR[k] < t1)
					s1-=(t1-ipR[k]);
				if (ipD[k] > t2)
					s2-=(ipD[k]-t2);
				d=(s1 < s2)? s1: s2;
				if (d > delta)
					d=delta;
				if (d > 0) {
					w+=(m_dpArray[sz]=d)*dpX[k];
					m_ipArray[sz++]=k;
					if (ipR[k] > m_ipR[k] || ipD[k] < m_ipD[k])
						type|=CTR_LOCAL;
				}
			} // for (int k=sz=0;
			if (w > 1.01*delta) {
				if (genFlag) {
					safeAddCut(-2,type,-INF,delta,sz,m_dpArray,m_ipArray);
//					printRow(m_iM-1);
					++cutNum;
					ipR=reinterpret_cast<int*>(m_dpFb); ipD=ipR+n;
				}
				else return true;
			}
		} // for (int j=0;
	} // for (int i=0;
	return (cutNum)? true: false;
} // end of COneMachine::separate()

/////////////////////////////////////////////////////
void COneMachine::printSchedule(const char* fileName)
{
    std::ofstream fout(fileName);
    if (!fout.is_open()) {
        throw new CFileException("COneMachine::printSchedule",fileName);
    }
    double *dpX;
    int n,s,*ipHd;
    getSolution(dpX,ipHd);
    n=m_iJobNum;
    fout << "task (c) [r,d;p]: start end\n";
    for (int i=0; i < n; i++) {
        if (dpX[i] > 0.5) {
            s=(int)dpX[n+i];
            fout << i << " (" << getObjCoeff(i) << ") "
                << " [" << m_ipR[i] << "," << m_ipD[i] << ";" << m_ipP[i] << "]: "
                << s << ", " << s+m_ipP[i] << std::endl;
        }
    }
    fout.close();
}
