#include <cmip.h>

class CCutStock: public CMIP
{
	double *m_dpF;
	int	*m_ipFinalLength;
	int *m_ipFinalNum;
	int m_iRawLength;
	int m_iFinalTypeNum;

public:
	CCutStock(const char* name);
	virtual ~CCutStock();

	#ifndef __ONE_THREAD_
	CCutStock(const CCutStock &other, int thread);
	CMIP* clone(const CMIP *pMip, int thread);
	#endif

	void init(int m, int *l, int *q, int L);
	void newColumn(int &k, int *S, int *b, int &sz);

	bool generateColumns(int m, const tagHANDLE *ipRowHd, const double *dpY) final;
	static double intKnapsack(int n, double *c, int *a, int b, int*x, double *dpMem=0);

	void printSolution(const char* name=0);
};
