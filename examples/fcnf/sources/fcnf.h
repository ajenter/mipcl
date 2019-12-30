#include <cmip.h>
#include <except.h>

class CFCNF: public CMIP
{
	bool m_bMemory;
	int m_iVertNum,m_iEdgeNum;
	int *m_ipHead, *m_ipTail,
		*m_ipFxCost, *m_ipCost,
		*m_ipCap,*m_ipDem;
public:
	CFCNF(const char* name): CMIP(name), m_bMemory(false), m_ipTail(0)
	{};  ///< constructor 1

	CFCNF(const char* name,
		int n, int m, int* tail, int* head,
		int *capacity, int* fixedCost, int* cost,
		int* demand): CMIP(name), m_bMemory(false),
			m_iVertNum(n), m_iEdgeNum(m),
			m_ipTail(tail), m_ipHead(head),
			m_ipCap(capacity),
			m_ipFxCost(fixedCost), m_ipCost(cost),
			m_ipDem(demand)
	{};  ///< constructor 2

	virtual ~CFCNF() {if (m_bMemory) delete[] m_ipTail;}; ///< destructor

	// implementation
	void readNet(const char* fileName);
	void model();
	void printSolution(const char* fileName); ///< overrides base class function
};
