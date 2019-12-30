// tsp.h: interface for the Ctsp class.
//////////////////////////////////////////////////
#ifndef __TSP__H
#define __TSP__H

#include <cmip.h>
#include "flowNet.h"

class CTspPool;
enum EDGE_WEIGHT_TYPE {ATT,GEO,EUC_2D};

class CTsp: public CMIP  
{
	CFlowNet<double> *m_pNet; ///< implements graph algorithms used for generasting cuts
	EDGE_WEIGHT_TYPE m_enType; ///< specifies how distances are computed
	int m_iPointNum;	///< number of points;
	double* m_dpCoordX; ///< x-coordinates
	double* m_dpCoordY; ///< y-coordinates
	double m_dMaxDist; ///< max distance between two point
	double m_dTourLength; ///< length of best tour found so far
	int *m_ipNextOnTour; ///< best tour
	CTspPool *m_pTspPool; ///< pool for storing cuts

public:
	CTsp(const char *name); ///< base constructor
	#ifndef __ONE_THREAD_
		CTsp(const CTsp &other, int thread);
		CMIP* clone(const CMIP *pMip, int thread);
	#endif
		
	virtual ~CTsp(); ///< destructor

private:
	double dist(int i, int j); ///< computes the distance between points `i` and `j`

	void allocMemForCoords();

	void readPoints(const char *fileName); ///< reads coordinates of points
	void buildActiveGraph();
	void setMIP();

    bool getRow(tagHANDLE hd, int n, const tagHANDLE* ipColHd,
                        int& type, double& b1, double& b2,
                        int& sz, double* dpVal, int* ipCol, bool &bScaled);
    bool getColumn(tagHANDLE hd, int m, const tagHANDLE* ipRowHd,
                int& type, double& cost, double& d1, double& d2,
                int& sz, double* dpVal, int* ipRow);
	void lockCtr(tagHANDLE hd);
	void unlockCtr(tagHANDLE hd);

	double goToNearest(int s, int* tau); // builds initial feasible solution
	void approximate(); // calls goToNearest() in turns with any node as a starting position

// separation routines
	bool separate(int n, const double* dpX, const tagHANDLE* ipColHd, bool GenFlag);
	int separateFromPool(int n, const double* dpX, const int* ipColHd, bool bGenFlag);
	void buildSupportGraph(int n, const double* dpX, const tagHANDLE* ipColHd);
	int cutSeparate(int varNum, const double* dpX, const tagHANDLE* ipColHd, bool GenFlag);

	bool genCut1(int n, const double* dpX, const tagHANDLE* ipColHd);
	int BLOSSOM_separate(int n, const double* dpX, const tagHANDLE* ipColHd);
	void addCombCut(int b, int* ipComb,
		int n, const double* dpX, const tagHANDLE* ipColHd);

	bool generateColumns(int m, const tagHANDLE* ipRowHd, const double* dpY);
	void buildSupportGraph();

	void changeRecord(double dObjVal, int n, const double* dpX, const tagHANDLE* ipHd);
public:
	void solve();
	void printSolution(const char* fileName=0);
}; //==============================================

#endif // #ifndef __TSP__H
