#include <except.h>
#include <cmip.h>

class Coptportfolio: public CMIP
{
	int m_iN, m_iT, m_iQ1, m_iQ2;
	double m_dV;
	double  *m_dpRet, *m_dpMu, *m_dpP, *m_dpL;
public:
	Coptportfolio(const char* name);
#ifndef __ONE_THREAD_
	Coptportfolio(const Coptportfolio &other, int thread);
	CMIP* clone(const CMIP *pMip, int thread);
#endif
	virtual ~Coptportfolio();
//////
	void readData(const char* fileName);
	void printSolution(const char* fileName) ;

//////
	void model(double p, double V);
private:
	void computeParameters(double p);
};

