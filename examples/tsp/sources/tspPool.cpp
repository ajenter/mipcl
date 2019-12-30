#include <iostream>
#include <cstring>
#include <new>
#include <except.h>
#include "tspPool.h"

CTspPool::CTspPool(int ptNum): m_iPointNum(ptNum) 
{
	_MUTEX_ASSIGN_P(m_pMemMutex,0)
	m_iMaxPoolSize=10*ptNum;
	m_iPoolSize=0;
	m_iFirstFreeInPool=-1;
	allocMem();
	_RWLOCK_INIT(m_rwLock)
}

void CTspPool::allocMem()
{
	m_pEntry = new(std::nothrow) tagPoolEntry[m_iMaxPoolSize];
	m_ipBuf=new(std::nothrow) int[m_iMaxPoolSize*m_iPointNum];
	if (!m_pEntry || !m_ipBuf ) {
		throw new CMemoryException("CTspPool::allocMem");
	}
} // end of CTspPool::allocMem

CTspPool::~CTspPool()
{
	if (m_pEntry)
		delete[] m_pEntry;
	if (m_ipBuf)
		delete[] m_ipBuf;
}
	
void CTspPool::reallocMem()
{
	int *ipBuf;
	tagPoolEntry *pEntry=0;
	int sz=m_iMaxPoolSize+m_iMaxPoolSize/3;
	_MUTEX_SAFE_LOCK(m_pMemMutex)
	pEntry=new(std::nothrow) tagPoolEntry[sz];
	ipBuf=new(std::nothrow) int[sz*m_iPointNum];
 	_MUTEX_SAFE_UNLOCK(m_pMemMutex)
	if (!ipBuf || !pEntry) {
		throw new CMemoryException("CTspPool::reallocMem");
	}
	memcpy(ipBuf,m_ipBuf,m_iMaxPoolSize*m_iPointNum*sizeof(int));
	for (int i=0; i < m_iMaxPoolSize; ++i) {
		pEntry[i].state=m_pEntry[i].state;
		pEntry[i].ct=m_pEntry[i].ct;
		pEntry[i].b=m_pEntry[i].b;
	}
	_MUTEX_SAFE_LOCK(m_pMemMutex)
	delete m_pEntry;
	delete m_ipBuf;
 	_MUTEX_SAFE_UNLOCK(m_pMemMutex)
	m_ipBuf=ipBuf;
	m_pEntry=pEntry;
	m_iMaxPoolSize=sz;
} // end of CTspPool::reallocMem

////////// P O O L  handling

int CTspPool::getCoefficient(int yv, int yw)
{
	int s1,s2,l=0;
	s1=yv & 0x0000FFFF;
	s2=yw & 0x0000FFFF;
	if (s1 && (s1==s2))
		l++;
	s1=yv >> 16;
	s2=yw >> 16;
	if (s1 && (s1==s2))
		l++;
	return l;
} // end of Ctsp::getCoefficient

int CTspPool::addCut(int b, int *ipComb)
{
	int hd;
	_RWLOCK_WRLOCK(&m_rwLock)
	if (m_iFirstFreeInPool < 0) {
		if (m_iPoolSize >= m_iMaxPoolSize) {
			freeNotUsedCuts();
			if (m_iFirstFreeInPool < 0)
				reallocMem();
		}
	}
	if (m_iFirstFreeInPool < 0)
		hd=m_iPoolSize++;		
	else  {
		hd=m_iFirstFreeInPool;
		m_iFirstFreeInPool=m_pEntry[hd].ct;
	}
	m_pEntry[hd].state=0;
	m_pEntry[hd].ct=1;
	m_pEntry[hd].b=b;
	memcpy(m_ipBuf+(hd*m_iPointNum),ipComb,m_iPointNum*sizeof(int));
	_RWLOCK_UNLOCK_WRLOCK(&m_rwLock)
	return hd;
} // end of CTspPool::addCut

void CTspPool::freeNotUsedCuts()
{
	tagPoolEntry *pEntry=m_pEntry;
	for (int i=0; i < m_iMaxPoolSize; ++i, ++pEntry) {
		if (!pEntry->state && !pEntry->ct) {
			pEntry->state=-1;
			pEntry->ct=m_iFirstFreeInPool;
			m_iFirstFreeInPool=i;
		}
	}
} // end of CTspPool::freeNotUsedCuts

int CTspPool::buildRow(int hd, int n, const int *ipColHd,
				   double* dpVal, int* ipCol, double &rhs)
{
	int l,colHd, sz, *data;
	_RWLOCK_RDLOCK(&m_rwLock)
	data=m_ipBuf+(hd*m_iPointNum);
	for (int i=sz=0; i < n; ++i) {
		colHd=ipColHd[i];
		if  (l=getCoefficient(data[colHd >> 16],data[colHd & 0x0000FFFF])) {
			dpVal[sz]=l;
			ipCol[sz++]=i;
		}
	}
	rhs=m_pEntry[hd].b;
	_RWLOCK_UNLOCK_RDLOCK(&m_rwLock)
	return sz;
} // end of CTspPool::buildRow

int CTspPool::buildColumn(int v, int w, int m, const int* ipRowHd,
					double* dpVal, int* ipRow)
{
	int l,*data,sz=2;
	dpVal[0]=dpVal[1]=1.0;
	ipRow[0]=v; ipRow[1]=w;
	_RWLOCK_RDLOCK(&m_rwLock)
	for (int i=m_iPointNum; i < m; ++i) {
		if ((l=ipRowHd[i]) >= 0) {
			data=m_ipBuf+(l*m_iPointNum);
			if (l=getCoefficient(data[v],data[w])) {
				dpVal[sz]=l;
				ipRow[sz++]=i;
			}
		}
	}	// for (register int i=0;
	_RWLOCK_UNLOCK_RDLOCK(&m_rwLock)
	return sz;
} // end of CTsp::buildColumn

int CTspPool::getNextCut(int thread, int &hd,
		int n, const double* dpX, const int* ipColHd,
		int *ipCol, double *dpVal, double &rhs)
{
	double w;
	int sz, colHd, a, *data;
	tagPoolEntry *pEntry;
	sz=0;
	_RWLOCK_RDLOCK(&m_rwLock)
	data=m_ipBuf+(hd*m_iPointNum);
	pEntry=m_pEntry+hd;
	for (; hd < m_iPoolSize; ++hd, ++pEntry) {
		if (pEntry->state >= 0 && !(pEntry->state & (1 << thread))) {
			w=pEntry->b;
			for (int e=0; e < n; ++e) {
				colHd=ipColHd[e];
				a=getCoefficient(data[colHd >> 16],data[colHd & 0x0000FFFF]);
				if (a > 0) {
					w-=dpX[e]*(dpVal[sz]=a);
					ipCol[sz++]=e;
				}
			}
			if (w < -0.001) {
				rhs=pEntry->b;
				break;
			}
			sz=0;
		}
		data+=m_iPointNum;
	}
	_RWLOCK_UNLOCK_RDLOCK(&m_rwLock)
	return sz;
} // end of CTspPool::getNextCut
