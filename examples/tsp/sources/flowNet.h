///////////////////////////////////////////////////////////////
/**
 * \file flowNet.h interface for CMIP class
 * |  __Author__  | N.N. Pisaruk                              |
 * |-------------:|:------------------------------------------|
 * | __Address__  | Belarus, 220020 Minsk, L.Ukrainky, 8, 166 |
 * |    __Tel.__  | (00 375 17)250-08-32                      |
 * |  __Mobile__  | +37529 2761930                            |
 * |  __e-mail__  | pisaruk@yandex.by                         |
 * | __home page__| http://pisaruk.narod.ru                   |
 *
 *   \copyright  __N.N. Pisaruk,      1995 - 2014__
 */
#ifndef __FLOWNET__H
#define __FLOWNET__H

enum tagProblem {
	undir_min_cut=0, gomory_tree=1, odd_cut=2,
	dir_max_flow=3, undir_max_flow=4, all_pairs_cuts=5,
	undir_transhipment=6, dir_transhipment = 7,
	feasible_solution=8, comb=9, bi_comp=10
};

template <typename Tflow> class CFlowNet
{
	bool m_bGrOK; ///< The flag is set to `true` if the initializer for a particular problem succeeded.
	void* m_pGrMemBuf; ///< Memory buffer for storing graphs.
	void* m_pMC_MemBuf; ///< Memory buffer for storing numeric parameters of minimum cut problem instances.
	void* m_pPR_MemBuf; ///< Memory buffer for storing numeric parameters of maximum flow problem instances.
	void* m_pGH_MemBuf;
	void* m_pTRAN_MemBuf; ///< Memory buffer for storing numeric parameters of transportation problem instances.

protected:
	int m_iVertNum; ///< number of vertices in the graph given by `m_ipTail` and `m_ipHead`.
	int m_iMaxVertNum; ///< maximum number of vertices in any graph that can be stored in this class without reallocating memory.
    int m_iEdgeNum; ///< number of edges in the graph given by `m_ipTail` and `m_ipHead`.
    int m_iMaxEdgeNum; ///< maximum number of edges in any graph that can be stored in this class without reallocating memory.
	int *m_ipTail; ///< array of size `m_iEdgeNum`, `m_ipTail[e]` is the tail of edge `e`.
	int *m_ipHead; ///< array of size `m_iEdgeNum`, `m_ipHead[e]` is the head of edge `e`
	int *m_ipEdge; ///< array of size `m_iEdgeNum`, list of edges.
	int *m_ipEdgeSep; ///< array of size `m_iVertNum+1`, where `m_ipEdge[m_ipEdgeSep[v]],...,m_ipEdge[m_ipEdgeSep[v+1]-1]` is set ov edges incident to vertex `v`.

	Tflow *m_ipUCap; ///< array of size `m_iEdgeNum`, `m_ipUCap[e]` is upper capacity of edge `e`.
	Tflow *m_ipLCap; ///< array of size `m_iEdgeNum`, `m_iLCap[e]` is lower capacity of edge `e`.

	Tflow *m_ipCost; ///< array of size `m_iEdgeNum`, `m_ipCost[e]` is cost of edge `e`.
	Tflow *m_ipDemand; ///< array of size `m_iVertNum`, `m_ipDemand[e]` is demand at vertex `v`.

	int m_iTermNum; ///< number of terminal vertices.
	int *m_ipTerminal; ///< array of size `m_iTermNum` storing terminals.

	tagProblem m_eProblem; ///< index of the problem that is being solved.
	int m_iProblemMsk;  // = m_ipMask[m_eProblem]
	static const int m_ipMask[];
public:
	static const int m_iSourceThinkMsk;
	static const int m_iTerminalMsk;
	static const int m_iDirectMsk;
	static const int m_iLoCapMsk;
	static const int m_iUpCapMsk;
	static const int m_iCostMsk;
	static const int m_iDemandMsk;

	Tflow m_iBigM;  ///< Big M value.

// internal data structures for different problems
private:
////////////////////////////////////////////////////
//             2-connected Components
//
	int *m_ipBI_MemBuf;	///< Memory buffer needed by the procedure computed 2-connected components.
	int *m_ipNum, *m_ipBack, *m_ipStack, *m_ipCompStack;
	int m_iCompNum;	///< number of 2-connected components.

	/**
	 * Arrays to store 2-connected components.
	 * The set of vertices of component `k` (`0 <= k <= m_iCompNum`) is
	 * {m_ipComp[m_ipCompSep[k]],...,m_ipComp[m_ipCompSep[k+1]-1]}.
	 */
	int *m_ipComp, *m_ipCompSep; ///< array to store 2-connected components

////////////////////////////////////////////////////
//           Minimum cuts in undirected graph
//
	Tflow*  m_ipCutCap;
	int*    m_ipName;
	int*    m_ipNext;
	bool*   m_bpFlag;
	int*   m_ipCut;
	Tflow   m_iCutVal;

////////////////////////////////////////////////////
//   Maximum flows (the Push and Relabel method
//
	Tflow* m_ipFlow; ///< array of size `m_iEdgeNum` to store arc flows.
	int m_iSource; ///< source vertex.
	int m_iThink; ///< think vertex.
	Tflow m_iMaxFlow; ///< value of maximum flow.

///< internal arrays
	int*   m_ipLabel;
	Tflow* m_ipExcess;
	int*   m_ipCurEdge;

///< Queue of size `m_iVertNum`
	int* m_ipQueue; ///< array of size `m_iVertNum` to store queue elements
	int  m_iFirst; ///< first element in queue
	int  m_iLast; ///< last element in queue

/////////////////////////////////////////////////////
//                   Gomory---Hu trees
// Implementation of the algorithm presented in
// "Dan Gusfield. Very simple methods for all pairs network flow analysis. SIAM J. Comput. Vol. 19, No.1, pp. 143-155
//
	int* m_ipP;  ///< for i=0,...,m_iVertNum-1, (i,p[i]) is in the Gomory---Hu cut tree
	Tflow* m_ipFl; ///< and has value of `ipFl[i]`
/////////////////////////////////////////////////////
//                  Minimum T-cut
//

/////////////////////////////////////////////////////
//               Transportation Problem
//
	Tflow* m_ipPrice; ///< array of size `m_iVertNum`
	int* m_ipPrec;
	int* m_ipDepth;
	int* m_ipThread;

	bool m_iTRAN_Feasible; ///< if set to `true`, the transportation problem has a feasible solution.

//    Common Functions
public:
	CFlowNet(); ///< The default constructor.


	CFlowNet(int maxVertNum, int maxEdgeNum, int problemMask);
	virtual ~CFlowNet();

	bool isOK()
		{return m_bGrOK;}

	///< \return number of vertices.
	int getVertNum()
		{return m_iVertNum;}

	///< \return maximum number of vertices in graph that can be stored in this object without reallocating memory.
	int getMaxVertNum()
		{return m_iMaxVertNum;}

	void setVertNum(int vertNum)
		{m_iVertNum=vertNum;}

	int getEdgeNum()
		{return m_iEdgeNum;}
	int getMaxEdgeNum()
		{return m_iMaxEdgeNum;}

	void setEdgeNum(int edgeNum)
		{m_iEdgeNum=edgeNum;}

	void setMask(int iMask)
		{m_iProblemMsk=iMask;}

	void setUpCap(Tflow* ipUpCap)
		{m_ipUCap=ipUpCap;}

	void setLoCap(Tflow* ipLoCap)
		{m_ipLCap=ipLoCap;}

	void setCost(Tflow* ipCost)
		{m_ipCost=ipCost;}

	void SsetDemand(Tflow* ipDemand)
		{m_ipDemand=ipDemand;}

	void reset(int n);	// n - number of vertices

	bool findFeasibleSol();

	void printFlow();
	void printCut();

	void setBigM(Tflow BigM)
		{m_iBigM=BigM;}

	tagProblem getProblem()
		{return m_eProblem;}

	void setProblem(tagProblem Probl)
		{m_iProblemMsk = m_ipMask[int(m_eProblem=Probl)];}


	void delLastEdge()
		{m_iEdgeNum--;}

	void delEdge(int Edg, bool bLast=false);

	void GrgetMem(int MaxVertNum=0, int MaxEdgeNum=0);
	bool addEdge(int v, int w);
	bool addEdge(int v, int w, Tflow cap);
	int  getEdgeNo(int v, int w);

	///< The procedure builds the list of incident edges.
	void buildEdgeList();

	void serialize(char *GrpName);
// read the input file

	void printNet();
////////////////////////////////////////////////////
//                 Disjoint Sets                  //
////////////////////////////////////////////////////
private:
	inline void DS_MakeSet(int x);
	int DS_Find(int x);
	int DS_Link(int x, int y);

////////////////////////////////////////////////////
//              Connected Components              //
////////////////////////////////////////////////////
public:
	///< The procedure allocates memory for running the procedure that computes 2-connected components.
	void BI_getMem();

	///< The procedure releases memory allocated by `BI_getMem()`.
	void BI_freeMem();

	///< The procedure that computes 2-connected components.
	void BI_Comp();

	///< \return number of 2-connected components.
	int BI_getCompNum()
		{return m_iCompNum;}

    /**
     * \param[in] i index of a component, `i=0,...,m_iCompNum-``.
     * \param[out] ipComp array to store vertices that belong to component `i`.
     * \return size of Component `i`.
     */
	int BI_getComp(int i, int*& ipComp);

	///< The procedure prints all 2-connected components.
	void BI_printSolution();

////////////////////////////////////////////////////
//      Minimum Cut in undirected graphs
//
	///< The procedure allocates memory for running the procedure that computes a minimum cut in undirected graphs.
	void MC_getMem();

	///< The procedure releases memory allocated by `MC_getMem()`.
	void MC_freeMem();

	///< The procedure computes minimum cuts in undirected graphs.
	void MC_MinCut(Tflow m_iThreshold=0);

	///< \return the value of the cut found by `MC_MinCut()`.
	Tflow MC_getCutValue()
		{return m_iCutVal;}

	///< \return pointer to the array storing a description of the minimum cut found by `MC_MinCut()`.
	int* MC_getCut()
		{return m_ipCut;}

	void MC_printSolution(); ///< Prints the minimum cut found by `MC_MinCut()`.

////////////////////////////////////////////////////
//                Push and Relabel
//
	///< The procedure allocates memory for running the procedure that computes maximum flow in digraphs.
	void PR_getMem();

	///< The procedure releases memory allocated by `PR_getMem()`.
	void PR_freeMem();

	///< The procedure computes maximum flows in digraphs.
	void PR_pushAndRelabel();

	int  PR_BDF(); ///< \return the cardinality of the base set of a minimum cut.
	void PR_printSolution(); ///< Prints a maximum flow computed by `PR_pushAndRelabel()`.

///////////////////////////////////////////////////
//             Gomory---Hu Cut Tree
//
	void GH_getMem();
	void GH_freeMem();
	void GH_GomoryTree();
	void GH_printSolution();

///////////////////////////////////////////////////
//                 Minimum T-Cut
//
	void TCUT_getMem();
	void TCUT_freeMem();
	void T_Cut();

///////////////////////////////////////////////////
//                 Minimum Combs
//
	///< The procedure allocates memory for running the procedures generating combs.
	void COMB_getMem();

	///< The procedure releases memory allocated by `COMB_getMem()`.
	void COMB_freeMem();

////////////////////////////////////////////////////////////////////////////////
// ATTENTION: This procedure is designed to be used in a certain environment, //
//            it is assumed that when allocating memory,                      //
//            m_iMaxVertNum >= m_iVertNum + m_iEdgeNum                        //
//            m_iMaxEdgeNum >= 2*m_iEdgeNum                                   //
//            Be very careful when using it!                                 //
////////////////////////////////////////////////////////////////////////////////
	void blossom();

/**
 *  The procedure is called iteratively to extract combs computed by `blossom()`.
 * \param[in,out] s iterator index, for the first call `s=1`;
 * \param[in] delta only those blossoms are returned which induces inequalities violated by more than `delta`;
 * \param[out] hdSize size of handle;
 * \param[out] ipComb pointer to array of size `m_iVertNum`, describes returned comb.
 * \return number of teeth, or 0 if there is no comb.
 */
	int getNextBlossom(int &s, double delta, int &hdSize, int* &ipComb);

	int maxBlossom(); ///< \return number of combs found by `blossom()`.

/**
 * \param[in] k comb to be returned;
 * \param[out] hdSize size of handle;
 * \param[out] thNum number of teeth;
 * \param[out] ipComb pointer to array of size `m_iVertNum`, describes returned comb.
 */
	void getMaxBlossom(int k, int &hdSize, int &thNum, int* &ipComb);

	void COMB_printSolution(); ///< The procedure prints combs computed by `blossom()`.

	void COMB_AND_CUT_getMem();

///////////////////////////////////////////////////
//             Network Simplex Method            //
///////////////////////////////////////////////////
// Input:  Network (`m_ipTail,m_ipHead,m_ipUCap,m_ipLCap,m_ipCost,m_ipDemand`)
// Output: either
//         1) Cut `(S,V-S)` violating Hoffman's criterion
//         or
//         2) optimal flow m_ipFlow and optimal prizes `m_ipPrize`


	/**
	 *  The procedure allocates memory for needed to solve the transportation problem
	 *  in the network with the parameters stored in the arrays (`m_ipTail,m_ipHead,m_ipUCap,m_ipLCap,m_ipCost,m_ipDemand`).
	 */
	void TRAN_getMem();

	///< The procedure releases memory allocated by `TRAN_getMem()`.
	void TRAN_freeMem();

	/**
	 * The procedure implements the network simplex method.
	 * \param[in] initSol if `true`, an initial feasible solution is tored in
	 */
	void TRAN_simplex(bool initSol=false);
	void TRAN_initSolution();
	void TRAN_printSolution();
	void TRAN_tree();

///////////////////////////////////////////////////
//             Implementation                    //
///////////////////////////////////////////////////
	/**
	 * \param[in] e edge index;
	 * \return tail of edge `e`.
	 */
	int getTail(int e)
		{return m_ipTail[e];}

	/**
	 * \param[in] e edge index;
	 * \return head of edge `e`.
	 */
	int getHead(int e)
		{return m_ipHead[e];}

	/**
	 * \param[in] v vertex index;
	 * \param[out] ipList list (array) of edges incident to `v`.
	 * \return number of edges storing in array `ipList`.
	 */
	int	getEdgeList(int v, int* ipList);

	/**
	 * \param[in] e edge index;
	 * \return lower capacity of edge `e`.
	 */
	Tflow getLoCap(int e);

	/**
	 * \param[in] e edge index;
	 * \return upper capacity of edge `e`.
	 */
	Tflow getEdgeUCap(int e)
		{return m_ipUCap[e];}

private:

	Tflow resCap(int e, int v);

	/**
	 * \param[in] e edge index;
	 * \param[in] v vertex incident to `e`.
	 * \return `w`, where `e=(v,w)` or `e=(w,v)`.
	 */
	int getOtherEnd(int e, int v);

///////////////////////////////////////////////////
//      Push and Relabel: Implementation
//
	void initQueue();
	int  getFromQueue();
	void addToQueue(int v);
	void PR_initSolution();
	void PR_discharge(int v);
	void PR_relabel(int v);
	void PR_push(int v, int e);
	int	 PR_getNextFeasible(int v, int &i);

//////////////////////////////////////////////////
//        Network Simplex: Implementation
//

	int  TRAN_checkFlowOpt(int& direct);
	int  TRAN_changeFlow(int edge, int& side, int direct);
	void TRAN_computePrices();
	void TRAN_initTree();
	void getEmptyTree();
	int  TRAN_findRoot(int vertex);
	void TRAN_addEdge(int edge, int direct=1);
	int  TRAN_delEdge(int edge, int delta=0);
	void TRAN_setBigM();

	void print();
};

#endif // #ifndef __FLOWNET__H
