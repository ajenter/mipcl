#include <cmip.h>

class CGenAssign: public CMIP
{
	bool m_bOK;
	int m_iTskNum, m_iMachNum; ///< number of tasks and number of machines
	int *m_ipMachCap, *m_ipCost, *m_ipProcTime; ///< only pointers
	int *m_ipAssign; ///< best assignment, task i is assign to machine m_ipAssign[i]
	int m_iSizeOfNodeData; ///< size of array m_ipNd
	/** m_ipNd is used for storing matrix $\Gamma$, where
	 * $\gamma_{ij}=0$, if task i is not assigned to any machine
	 * $\gamma_{ij}=1$, if task i is assigned to machine j
	 * $\gamma_{ij}=2$, if task i is not allowed to be assigned to machine j
	 */
	int* m_ipNd;
	/**
	 *  both are set in StartBranching
	 * the current node is divided into 2 branches:
	 * 1) task m_iCurTsk is assigned to machine m_iCurMach
	 * 2) task m_iCurTsk must not be assigned to machine m_iCurMach
	 */
	int m_iCurTsk, m_iCurMach;
	int m_iQ; ///< max machine capacity
	double *m_dpKn;
	int *m_ipKn;

public:
	CGenAssign(const char* name, int m, int n, int* l, int* c, int* p);
#ifndef __ONE_THREAD_
	CGenAssign(const CGenAssign &other, int thread);
	CMIP* clone(const CMIP *pMip, int thread);
#endif	
	virtual ~CGenAssign();
private:
	void buildMaster(); ///< builds master problem
	void setGammaEntry(int i, int j, int val);
	int getGammaEntry(int i, int j);
	bool generateColumns(int ctrNum, const tagHANDLE* ipRowHd, const double* dpY);
	int startBranching(int /*nodeHeight*/);
	bool updateBranch(int k);
	int storeNodeData(void*& pMem);
	void restoreNodeData(void* pMem);
	void changeRecord(double objVal, int varNum, const double* dpX, const tagHANDLE* ipHd);
public:
	void printSolution(const char* fileName);
};
