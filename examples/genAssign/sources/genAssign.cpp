#include <fstream>
#include <cmath>
#include <cstring>
#include <except.h>
#include "knapsack.h"
#include "genAssign.h"

CGenAssign::CGenAssign(const char* name, int m, int n, int* l, int* c, int* p):
	CMIP(name), m_bOK{true},
	m_iTskNum(m), m_iMachNum(n),
	m_ipMachCap(l), m_ipCost(c), m_ipProcTime(p)
{
	int q{m*n};
	int k{q>>4};
	if (q & 0xF)
		++k;
	q=0;
	for (int j={0}; j < n; j++) {
        if (l[j] > q)
            q=l[j];
    }
	try {
		m_ipAssign= new(std::nothrow) int[m];
		m_ipNd = new(std::nothrow) int[m_iSizeOfNodeData=k];
		KNAPSACK::getMemForBinKnapsack(m,m_iQ=q,m_dpKn,m_ipKn);
	} catch(std::bad_alloc* e) {
		m_bOK=false;
		strcpy(m_sWarningMsg,"CGenAssign::CGenAssign(: ");
		strncat(m_sWarningMsg,e->what(),128);
		delete e;
    }
	if (m_bOK) {
		memset(m_ipNd,0,m_iSizeOfNodeData*sizeof(int));
		buildMaster();
	}
} // end of CGenAssign::CGenAssign()

#ifndef __ONE_THREAD_
CGenAssign::CGenAssign(const CGenAssign &other, int thread): CMIP(other,thread)
{
	m_bOK=true;
	m_iTskNum=other.m_iTskNum;
	m_iMachNum=other.m_iMachNum;
	m_ipMachCap=other.m_ipMachCap;
	m_iQ=other.m_iQ;
	m_ipCost=other.m_ipCost;
	m_ipProcTime=other.m_ipProcTime;
	m_ipAssign=other.m_ipAssign;
	try {
		m_ipNd=new(std::nothrow) int[m_iSizeOfNodeData=other.m_iSizeOfNodeData];
		KNAPSACK::getMemForBinKnapsack(m_iTskNum,m_iQ,m_dpKn,m_ipKn);
	} catch(std::bad_alloc* e) {
		m_bOK=false;
		strcpy(m_sWarningMsg,"CGenAssign::CGenAssign(: ");
		strncat(m_sWarningMsg,e->what(),128);
		delete e;
    }
	if (m_bOK)
		memset(m_ipNd,0,m_iSizeOfNodeData*sizeof(int));
}

CMIP* CGenAssign::clone(const CMIP *pMip, int thread)
{
	return static_cast<CMIP*>(new CGenAssign(*static_cast<CGenAssign*>(const_cast<CMIP*>(pMip)),thread));
}
#endif

CGenAssign::~CGenAssign()
{
#ifndef __ONE_THREAD_
	if (!m_iThread) {
#endif
    if (m_ipAssign) {
        delete[] m_ipAssign;
        m_ipAssign=0;
    }
#ifndef __ONE_THREAD_
	}
#endif
	if (m_ipNd) {
        delete[] m_ipNd;
        m_ipNd=0;
	}
	if (m_dpKn) {
		delete[] m_dpKn;
		m_dpKn=0;
	}
	if (m_ipKn) {
		delete[] m_ipKn;
		m_ipKn=0;
	}
}

void CGenAssign::buildMaster()
{
	int *l{m_ipMachCap}, *c{m_ipCost};
	int bigM{0}, m{m_iTskNum}, n{m_iMachNum};
	for (int i{0}; i < m; i++) {
		int q{0};
		for (int j{0}; j < n; ++j) {
			if (q < c[j*m+i])
				q=c[j*m+i];
		}
		bigM+=q;
	}
    double one{1.0};
    openMatrix(m+n,m,m,false,true,0,5*m,5*n*m);

    for (int i{0}; i < m; ++i) {
        addCtr(i,0,1.0,1.0);
    }

    for (int i{0}; i < n; ++i) {
        addCtr(m+i,0,-INF,1.0);
    }

    for (int i{0}; i < m; ++i) {
        addColumn(i,0,-bigM,0.0,1.0,1,&one,&i);
    }
    preprocOff();
    setScaling(SCL_NO);
    closeMatrix();
} // end of GenAssign::buildMaster()

void CGenAssign::setGammaEntry(int i, int j, int val)
{
	int k0{i*m_iTskNum};
	int k{k0+j};
	m_ipNd[k >> 4]|= (val << ((k & 0xF)<<1));
	if (val == 1) {
		for (int s{0}; s < m_iMachNum; ++s)
			if (s != j) {
				k=k0+s;
				m_ipNd[k>>4]|= (2 << ((k & 0xF)<<1));
			}
	}
} // end of CGenAssign::setGammaEntry()

int CGenAssign::getGammaEntry(int i, int j)
{
	int k{j*m_iMachNum+i}, s{2*((k & 0xF)<<1)};
	return (m_ipNd[k>>4] & (3 << s)) >> s;
} // end of CGenAssign::getGammaEntry()

int CGenAssign::storeNodeData(void*& pMem)
{
    pMem=m_ipNd;
    return m_iSizeOfNodeData*sizeof(int);
} // end of CGenAssign::storeNodeData()

void CGenAssign::restoreNodeData(void* pMem)
{
    int *ipColToMach{m_ipArray},
    	*ipCol{reinterpret_cast<int*>(m_dpArray)};
    int col,sz, m{m_iTskNum}, n{m_iMachNum};

    memcpy(m_ipNd,pMem,m_iSizeOfNodeData*sizeof(int));

    for (int j{0}; j < n; j++) {
        sz=getRow(j+m,0,ipCol);
        for (int i=0; i < sz; ++i)
            ipColToMach[ipCol[i]]=j;
    }

    for (int i{0}; i < m; ++i) {
		sz=getRow(i,0,ipCol);
		for (int j{0}; j < sz; ++j) {
			if ((col=ipCol[j]) >= m)
				if (getGammaEntry(i,ipColToMach[col]) == 2)
					setVarUpBound(col,0.0);
        }
    }
} // end of CGenAssign::restoreNodeData()

int CGenAssign::startBranching(int nodeHeight)
{
    double w,v,val, e{0.9999999};
    int sz, m{m_iTskNum}, n{m_iMachNum}, mach{-1},
        *ipCol=m_ipArray;
    m_iCurMach=m_iCurTsk=-1;
    for (int j{0}; j < n; j++) {
        sz=getRow(j+m,0,ipCol);
        v=w=0.0;
        for (int i=0; i < sz; i++) {
            val=getVarValue(ipCol[i]);
            v+=val;
            w+=(val*val);
        }
        v=1.0-v;
        w+=(v*v);
        if (e > w) {
            e=w;
            mach=j;
        }
    }
    if (mach >= 0) {
        double *dpVal{m_dpArray};
        int sz1, *ipRow{reinterpret_cast<int*>(dpVal+m)};
        memset(dpVal,0,m*sizeof(double));
        sz=getRow(mach+m,0,ipCol);
        for (int i{0}; i < sz; ++i) {
            sz1=getColumn(ipCol[i],0,ipRow);
            val=getVarValue(ipCol[i]);
            for (int j{0}; j < sz1; ++j)
				if (ipRow[j] < m)
					dpVal[ipRow[j]]+=val;
        }
        w=2.0; sz=-1;
        for (int i{0}; i < m; ++i) {
            if (w > fabs(dpVal[i]-0.5)) {
                w=fabs(dpVal[i]-0.5);
                sz=i;
            }
        }
        m_iCurMach=mach;
        m_iCurTsk=sz;       
    }
    return (m_iCurMach >= 0 && m_iCurTsk >= 0)? 2: 0;
} // end of CGenAssign::startBranching()

bool CGenAssign::updateBranch(int k)
{
    int *ipCol{m_ipArray};
	bool *bpColToMach{reinterpret_cast<bool*>(m_dpArray)};
    int sz, m{m_iTskNum}, n{m_iMachNum};

    setGammaEntry(m_iCurTsk,m_iCurMach,k+1);
	memset(bpColToMach,0,getVarNum()*sizeof(bool));
    sz=getRow(m_iCurMach+m,0,ipCol);
	for (int j{0}; j < sz; ++j) {
        bpColToMach[ipCol[j]]=true;
    }
	sz=getRow(m_iCurTsk,0,ipCol);
	for (int j{0}; j < sz; ++j) {
		if ((k && !bpColToMach[ipCol[j]]) || (!k && bpColToMach[ipCol[j]]))
			setVarUpBound(ipCol[j],0.0);
    }
	return true;
} // end of CGenAssign::updateBranch()

bool CGenAssign::generateColumns(int ctrNum, const tagHANDLE* ipRowHd, const double* dpY)
{
    double c0, z, tol=getRedCostTol();
    int b, num, m{m_iTskNum}, n{m_iMachNum};
    double *w{m_dpArray+(m+1)};
	int *x, *a, *ipTsk, *ipRow{m_ipArray},
        *c{m_ipCost}, *p{m_ipProcTime}, *l{m_ipMachCap};
   ipTsk=ipRow +(m+1); x=(a=ipTsk+m)+m;

    for (int i{0}; i <= m; ++i) {
        m_dpArray[i]=1.0;
    }

    num=0;
    for (int j{0}; j < n; ++j) {
    	int q{0}, sz{0}, k{0};
        b=l[j];
        c0=0.0;
        for (int e, i{0}; i < m; ++i) {
            if (!(e=getGammaEntry(i,j))) {
				if (p[i] <= l[j]) {
					z=-c[i]-dpY[i];
					if (z > tol) {
						ipTsk[k]=i;
						w[k]=z;
						a[k++]=p[i];
					}
				}
            }
            else if (e == 1) {
                q+=c[i];
                ipRow[sz++]=i;
                b-=p[i];
                c0-=(c[i]+dpY[i]);
            }
        }
        if (k && b > 0) {
	        if (KNAPSACK::binKnapsack(k,w,a,b,x,m_dpKn,m_ipKn) > dpY[m+j]-c0+tol) {
				for (int i{0}; i < k; ++i) {
	                if (x[i]) {
	                    q+=c[ipTsk[i]];
	                    ipRow[sz++]=ipTsk[i];
	                }
	            }
	            ipRow[sz++]=m+j;
	            sz=addNewColumn(-1,VAR_BIN,-q,0.0,1.0,sz,m_dpArray,m_ipArray,false,false,0,true);
	            ++num;
			    w=m_dpArray+(m+1); ///< it is possible that memory for `m_ipArray` and `m_dpArray` has been reallocated
				x=(a=(ipTsk=(ipRow=m_ipArray)+(m+1))+m)+m;
			}
        }
        p+=m; c+=m;
    }
    return (num)? true: false;
} // end of CGenAssign::generateColumns()

void CGenAssign::changeRecord(double objVal,
			int varNum, const double* dpX, const tagHANDLE* ipHd)
{
	int *ipRow{m_ipArray};
    int m{m_iTskNum}, m0{m_iTskNum+m_iMachNum};
    for (int k,sz, j{m}; j < varNum; ++j) {
        if (dpX[j] > 0.5) {
            sz=getColumn(j,0,ipRow);
            for (int i{0}; i < sz; ++i) {
				if (ipRow[i] < m0 && ipRow[i] >= m) {
					k=ipRow[i]-m;
					break;
				}
			}
            for (int i{0}; i < sz; ++i) {
				if (ipRow[i] < m)
					m_ipAssign[ipRow[i]]=k;
            }
        }
    }
} // end of CGenAssign::changeRecord()

void CGenAssign::printSolution(const char* fileName)
{
	char *str{CLP::m_sWarningMsg};
	if (fileName)
		strcpy(str,fileName);
	else getProblemName(str);
	strcat(str,".sol");
    std::ofstream fout(str);
    if (!fout.is_open()) {
        throw new CFileException("CGenAssign::printSolution",str);
    }
	fout << "Obj=" << -(int)floor(getObjVal()+0.5) << std::endl;
    for (int i{0}; i < m_iTskNum; i++) {
		fout << i << " => " << m_ipAssign[i] << std::endl; 
    }
	fout.close();
} // end of CGenAssign::printSolution()
