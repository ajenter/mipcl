#include <cmip.h>

class COneMachine: public CMIP
{
    int m_iJobNum; ///< number of jobs
    /**
     * All three are arrays of size `m_iJobNum`,
     * resp., processing times, release and due dates of jobs.
     */
    int *m_ipP, *m_ipR, *m_ipD;
    double *m_dpCost; ///< job costs
    int m_iDmax; ///< maximum of release dates
    int m_iJob1, m_iJob2; ///< are set in `startBranching()` to be used in `updateBranch()`
public:
    COneMachine(const char* name,int n, double* c, int *p, int *r, int* d);
#ifndef __ONE_THREAD_
	COneMachine(const COneMachine &other, int thread);
	CMIP* clone(const CMIP *pMip, int thread);
#endif
	~COneMachine() {}
private:
    bool isFeasible(int varNum, const double* dpX, const tagHANDLE* ipColHd);
    int startBranching(int nodeHeight);
    bool updateBranch(int k);
    bool getRow(tagHANDLE hd, int varNum, const tagHANDLE* ipColHd,
                int& type, double& b1, double& b2,
                int& sz, double* dpVal, int* ipCol, bool &scaled);

    bool separate(int varNum, const double* dpX, const tagHANDLE* ColHd, bool GenFlag);

public:
    void printSchedule(const char* fileName);
};
