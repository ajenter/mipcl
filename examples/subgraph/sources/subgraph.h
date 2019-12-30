#include <cmip.h>
#include <except.h>

class CSubGraph: public CMIP
{
	bool m_bMemory;
	int m_iVertNum, m_iEdgeNum, m_iK;
	int *m_ipHead, *m_ipTail, *m_ipCost, *m_ipDeg;
	CLP* cut;
public:
	CSubGraph(const char* name): CMIP(name), m_bMemory(false), m_ipTail(0)
	{}
	CSubGraph(const char* name, int n, int m, int k,
			  int* ipTail, int* ipHead, int* ipCost, int* ipDeg): CMIP(name),
			  cut(0), m_iVertNum(n), m_iEdgeNum(m), m_iK(k), m_ipTail(ipTail), m_ipHead(ipHead),
			  m_bMemory(false)
	{}

// Clone constructor needed for multi-threading
#ifndef __ONE_THREAD_
	CSubGraph(const CSubGraph &other, int thread);
	CMIP* clone(const CMIP *pMip, int thread);
#endif

	virtual ~CSubGraph(void);

// implementation
	void readData(const char* fileName);
	void buildMatrix();
private:
	void cutInit(const char* name); // name is needed in order not to mix log streams
	bool separate(int varNum, const double* x, const tagHANDLE* ipColHd, bool bGenFlag);
};
