#ifndef TSPPOOL_H_
#define TSPPOOL_H_

///////////////////////////////////////////////////////////////////
// To store the cut lx <= l_0 having a left hand side of the form
// \sum_{i=1}^h x(E(H_i)) + \sum_{j=1}^t x(E(T_j)),
// H_i \cap H_k = \emptyset (i\neq k), T_j \cap T_r = \emptyset (j\neq r),
// we use a vector y  of length n=|V| defined as follows:
// each word y[v] is divided into 2 halfwords y_h[v] and y_l[v] where
// y_l[v] = i if  v\in H_i; and y_l[v] = 0 otherwise
// y_h[v] = j if  v\in T_j; and y_l[v] = 0 otherwise


#include <thread.h>

struct tagPoolEntry {
    int state; // number of node LPs that use this constraint
	int ct; // number of active nodes using this cut
	int b;
};

class CTspPool {
	int m_iPointNum, m_iPoolSize, m_iMaxPoolSize, m_iFirstFreeInPool; // describe the state of the pool
	tagPoolEntry* m_pEntry; // pool for storing cuts
	int* m_ipBuf; // memory buffer for storing cut coefficients
#ifndef __ONE_THREAD_
	_RWLOCK m_rwLock;
	_MUTEX *m_pMemMutex;
#endif
public:
	CTspPool(int ptNum);
	~CTspPool();
	void allocMem();
	void reallocMem();

	static int getCoefficient(int yv, int yw);
	int addCut(int b, int *ipComb);
	void freeNotUsedCuts();

	int buildRow(int hd, int n, const int *ipColHd, double* dpVal, int* ipCol, double &rhs);
	int buildColumn(int v, int w, int m, const int* ipRowHd, double* dpVal, int* ipRow);

	int getNextCut(int thread, int &hd, int n, const double* dpX, const int* ipColHd,
			int *ipCol, double *dpVal, double &rhs);

	void markCtr(int hd, int thread)
		{m_pEntry[hd].state|= 1 << thread;}
	void unmarkCtr(int hd, int thread)
		{m_pEntry[hd].state&=~(1 << thread);}
	
	void lockCtr(int hd)
	{
		_RWLOCK_WRLOCK(&m_rwLock)
		++m_pEntry[hd].ct;
		_RWLOCK_UNLOCK_WRLOCK(&m_rwLock)
	}
	void unlockCtr(int hd)
	{
		_RWLOCK_WRLOCK(&m_rwLock)
		--m_pEntry[hd].ct;
		_RWLOCK_UNLOCK_WRLOCK(&m_rwLock)
	}
	
	void wrLockPool()
		{_RWLOCK_WRLOCK(&m_rwLock)}
		
	void wrUnlockPool()
		{_RWLOCK_UNLOCK_WRLOCK(&m_rwLock)}

	void rdUnlockPool()
		{_RWLOCK_UNLOCK_RDLOCK(&m_rwLock)}

	int* getBufPtr()
		{return m_ipBuf;}

}; //========================================

#endif /*TSPPOOL_H_*/
