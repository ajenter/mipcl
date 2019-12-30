///////////////////////////////////////////////////////
// FlowNet.cpp: implementation of the FlowNet class. //
///////////////////////////////////////////////////////

#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
//#include <cstdio>
#include <except.h>
#include "compare.h"
#include "flowNet.h"

using std::cout;
using std::cout;
using std::endl;
using std::setw;

#define MAX_STR_LEN 256

const int NIL=-1;
const int DEL=-2;

template <class Tflow>
const int CFlowNet<Tflow>::m_ipMask[] = 
{0x0004,0x0004,0x0024,0x0015,0x0014,0x0004,0x0206,0x0307,0x0208,0x0004,0x0000};

template <class Tflow>
const int CFlowNet<Tflow>::m_iSourceThinkMsk = 0x0010;

template <class Tflow>
const int CFlowNet<Tflow>::m_iTerminalMsk = 0x0020;

template <class Tflow>
const int CFlowNet<Tflow>::m_iDirectMsk = 0x0001;

template <class Tflow>
const int CFlowNet<Tflow>::m_iLoCapMsk = 0x0002;

template <class Tflow>
const int CFlowNet<Tflow>::m_iUpCapMsk = 0x0004;


template <class Tflow>
const int CFlowNet<Tflow>::m_iCostMsk = 0x0100;

template <class Tflow>
const int CFlowNet<Tflow>::m_iDemandMsk = 0x0200;

//////////////////////////////////////////////////////////////////////
//               Construction/Destruction                           //
//////////////////////////////////////////////////////////////////////

template <class Tflow>
CFlowNet<Tflow>::CFlowNet()
{
	m_bGrOK=true;
	m_pGrMemBuf=0;
} // end of CFlowNet<Tflow>::CFlowNet

template <class Tflow>
CFlowNet<Tflow>::CFlowNet(int maxVertNum, int maxEdgeNum, int problemMask)
{
	m_iProblemMsk=problemMask;	
	GrgetMem(maxVertNum,maxEdgeNum);
}

template <class Tflow>
void CFlowNet<Tflow>::GrgetMem(int maxVertNum, int maxEdgeNum)
{
	if (maxVertNum)
		m_iMaxVertNum=maxVertNum;
	if (maxEdgeNum)
		m_iMaxEdgeNum=maxEdgeNum;
	int iSize=(2*sizeof(int)+sizeof(Tflow))*m_iMaxEdgeNum;
	if (m_iProblemMsk & m_iUpCapMsk)
		iSize+=(sizeof(Tflow)*m_iMaxEdgeNum);
	if (m_iProblemMsk & m_iLoCapMsk)
		iSize+=(sizeof(Tflow)*m_iMaxEdgeNum);
	if (m_iProblemMsk & m_iCostMsk)
		iSize+=(sizeof(Tflow)*m_iMaxEdgeNum);
	if (m_iProblemMsk & m_iDemandMsk)
		iSize+=(sizeof(Tflow)*m_iMaxVertNum);
	if (m_iProblemMsk & m_iTerminalMsk)
		iSize+=(sizeof(int)*m_iMaxVertNum);

	if (!(m_pGrMemBuf = new char[iSize])) {
		m_bGrOK=false;
		throw new CMemoryException("CFlowNet<Tflow>::GrgetMem");
	}

	m_ipHead =(m_ipTail=(int*)m_pGrMemBuf)+m_iMaxEdgeNum;
	Tflow* p=(Tflow*)(m_ipHead+m_iMaxEdgeNum);
	if (m_iProblemMsk & m_iUpCapMsk)
		p=(m_ipUCap=p)+m_iMaxEdgeNum;
	if (m_iProblemMsk & m_iLoCapMsk)
		p=(m_ipLCap=p)+m_iMaxEdgeNum;
	if (m_iProblemMsk & m_iCostMsk)
		p=(m_ipCost=p)+m_iMaxEdgeNum;
	if (m_iProblemMsk & m_iDemandMsk)
		p=(m_ipDemand = p)+m_iMaxVertNum;
	if (m_iProblemMsk & m_iTerminalMsk)
		m_ipTerminal=(int*)p;
}  // end of FlowNet::GrgetMem

template <class Tflow>
CFlowNet<Tflow>::~CFlowNet()
{
	if (m_pGrMemBuf)
		delete[] (Tflow*)m_pGrMemBuf;
}  // end of CFlowNet<Tflow>::~CGraph()

template <class Tflow>
void CFlowNet<Tflow>::serialize(char *GrpName)
{
	char fileName[MAX_STR_LEN];
	strcpy(fileName,GrpName);
	strcat(fileName,".grp");
	std::ifstream in_stream(fileName);
	if (!in_stream.is_open()) {
		m_bGrOK=false;
		throw new CFileException("CFlowNet<Tflow>::Serialize",fileName);
	}
	in_stream >> m_iVertNum;
	in_stream >> m_iEdgeNum;
	{
		int k; in_stream >> k; 
		m_eProblem=(tagProblem)k; 
		m_iProblemMsk=m_ipMask[k];
	} 

	m_iMaxVertNum= m_iVertNum;
	m_iMaxEdgeNum=m_iEdgeNum; 

	switch (m_eProblem) {
	//  case undir_max_flow:
	case feasible_solution:
	case dir_max_flow:
		if (m_iProblemMsk & m_iLoCapMsk)
			m_iMaxVertNum+=2;
		m_iMaxEdgeNum+=m_iVertNum; 
		break;
		//  case all_pairs_cuts:
		//      m_iMaxVertNum= m_iVertNum+2;
		//      m_iMaxEdgeNum=m_iEdgeNum+m_iVertNum;
		//	break;
	case undir_transhipment:
	case dir_transhipment:
		m_iMaxVertNum+=2; // enough +1 for Big M strategy
		m_iMaxEdgeNum+=m_iVertNum;
		break;
	case comb:
		m_iMaxVertNum+=m_iMaxEdgeNum;
		m_iMaxEdgeNum*=2;
		break;
	}

	GrgetMem(); // in case of failure throws CMemoryException

	for (int v=0; v < m_iEdgeNum; ++v) {
		in_stream >> m_ipTail[v];
		in_stream >> m_ipHead[v];

		// read problem specific data
		if (m_iProblemMsk & m_iLoCapMsk)
			in_stream >> m_ipLCap[v];
		if (m_iProblemMsk & m_iUpCapMsk)
			in_stream >> m_ipUCap[v];
		if (m_iProblemMsk & m_iCostMsk)
			in_stream >> m_ipCost[v];
	}	// end of for (register int v=0;

	if (m_iProblemMsk & m_iDemandMsk)
		for (int v=0; v < m_iVertNum;)
			in_stream >> m_ipDemand[v++];

	if (m_iProblemMsk & m_iSourceThinkMsk) {
		in_stream >> m_iSource;
		in_stream >> m_iThink;
	}

	if (m_iProblemMsk & m_iTerminalMsk) {
		in_stream >> m_iTermNum;
		for (int v=0; v < m_iTermNum; ++v)
			in_stream >> m_ipTerminal[v];
	}

	in_stream.close();
}  // end of CFlowNet<Tflow>::serialize

template <class Tflow>
void CFlowNet<Tflow>::reset(int n)
{
	m_iVertNum=n;
	m_iEdgeNum=0;
}	// end of CFlowNet<Tflow>::reset

template <class Tflow>
void CFlowNet<Tflow>::buildEdgeList()
{
	int m2=(m_iEdgeNum << 1);

	for (int v=0; v <= m_iVertNum; ++v)
		m_ipEdgeSep[v]=0;

	for (int e=0; e < m_iEdgeNum; ++e)
		if (m_ipHead[e] != DEL) {
			++m_ipEdgeSep[m_ipTail[e]+1];
			++m_ipEdgeSep[m_ipHead[e]+1];
		}

	m_ipEdgeSep[m_iVertNum]=m2-m_ipEdgeSep[m_iVertNum];

	for (int v=m_iVertNum-1; v >= 1; --v) {
		m_ipEdgeSep[v]=m_ipEdgeSep[v+1]-m_ipEdgeSep[v];
	}

	for (int e=0; e < m_iEdgeNum; ++e)
		if (m_ipHead[e] != DEL) {
			m_ipEdge[m_ipEdgeSep[m_ipHead[e]+1]++]=e;
			m_ipEdge[m_ipEdgeSep[m_ipTail[e]+1]++]=e;
		}
}  // end of CFlowNet<Tflow>::buildEdgeList

template <class Tflow>
int CFlowNet<Tflow>::getOtherEnd(int e, int v)
{
	return (m_ipTail[e] == v)? m_ipHead[e]: m_ipTail[e];
}  // end of CFlowNet<Tflow>::getOtherEnd

template <class Tflow>
int CFlowNet<Tflow>::getEdgeNo(int v, int w)
{
	int e;
	if (m_ipEdgeSep[v+1]-m_ipEdgeSep[v] > m_ipEdgeSep[w+1]-m_ipEdgeSep[w]) {
		e=v; v=w; w=e;
	}
	for (int i=m_ipEdgeSep[v]; i < m_ipEdgeSep[v+1]; ++i)
		if (m_ipTail[e=m_ipEdge[i]] != DEL)
			if (w==((m_ipTail[e] == v)? m_ipHead[e]: m_ipTail[e]))
				return e;
	return NIL;
}  // end of CFlowNet<Tflow>::getEdgeNo

template <class Tflow>
bool CFlowNet<Tflow>::addEdge(int v, int w)
{
	if (m_iEdgeNum == m_iMaxEdgeNum)
		return false;
	m_ipTail[m_iEdgeNum]=v;
	m_ipHead[m_iEdgeNum++]=w;
	return true;
}  // end of CFlowNet<Tflow>::addEdge

template <class Tflow>
bool CFlowNet<Tflow>::addEdge(int v, int w, Tflow cap)
{
	if (m_iEdgeNum == m_iMaxEdgeNum)
		return false;
	m_ipUCap[m_iEdgeNum]=cap;
	m_ipTail[m_iEdgeNum]=v;
	m_ipHead[m_iEdgeNum++]=w;
	return true;
}  // end of CFlowNet<Tflow>::addEdge

template <class Tflow>
void CFlowNet<Tflow>::delEdge(int Edg, bool bLast)
{
	m_ipTail[Edg]=DEL;
	if (bLast) {
		int edgeNo=0;
		for (int e=0; e < m_iEdgeNum; ++e) {
			if (m_ipTail[e] != DEL) {
				m_ipTail[edgeNo]=m_ipTail[e];
				m_ipHead[edgeNo]=m_ipHead[e];
				if (m_iProblemMsk & m_iLoCapMsk)
					m_ipLCap[edgeNo]=m_ipLCap[e];
				if (m_iProblemMsk & m_iUpCapMsk)
					m_ipUCap[edgeNo]=m_ipUCap[e];
				if (m_iProblemMsk & m_iCostMsk)
					m_ipCost[edgeNo]=m_ipCost[e];
				if (m_iProblemMsk & m_iDemandMsk)
					m_ipDemand[edgeNo]=m_ipDemand[e];
				++edgeNo;
			}	// if (m_ipTail[e] != DEL)
		}
		m_iEdgeNum=edgeNo;
	}	// if (bLast)
}  // end of CFlowNet<Tflow>::delEdge

template <class Tflow>
int CFlowNet<Tflow>::getEdgeList(int v, int* ipList)
{
	int k=m_ipEdgeSep[v+1]-m_ipEdgeSep[v];
	memcpy(ipList,m_ipEdge+m_ipEdgeSep[v],k*sizeof(int));
	return k;
}	// end of CFlowNet<Tflow>::getEdgeList

//////////////////////////////////////////////////////////////
//                 Connected Components                     //
//////////////////////////////////////////////////////////////
template <class Tflow>
void CFlowNet<Tflow>::BI_getMem()
{
	if (!(m_ipBI_MemBuf = new int[4*m_iMaxEdgeNum+6*m_iMaxVertNum+2])) {
		throw new CMemoryException("CFlowNet<Tflow>::BI_getMem");
	}
	m_ipEdgeSep=(m_ipEdge=m_ipBI_MemBuf)+(2*m_iMaxEdgeNum);
	m_ipNum=m_ipEdgeSep+(m_iMaxVertNum+1);
	m_ipBack=m_ipNum+m_iMaxVertNum;
	m_ipCurEdge=m_ipBack+m_iMaxVertNum;
	m_ipCompStack=(m_ipStack=m_ipCurEdge+m_iMaxVertNum)+m_iMaxVertNum;	
	m_ipCompSep=(m_ipComp=m_ipCompStack+m_iMaxEdgeNum)+m_iMaxEdgeNum;
}  // end of CFlowNet<Tflow>::BI_getMem

template <class Tflow>
void CFlowNet<Tflow>::BI_freeMem()
{
	if (m_ipBI_MemBuf)
		delete[] m_ipBI_MemBuf;
}  // end of  CFlowNet<Tflow>::BI_freeMem

template <class Tflow>
int CFlowNet<Tflow>::BI_getComp(int i, int*& Comp)
{
	if (i < 0 || i >= m_iCompNum)
		return NIL;
	Comp=m_ipComp+m_ipCompSep[i];
	return
			m_ipCompSep[i+1]-m_ipCompSep[i];
}	// end of CFlowNet<Tflow>::BI_getComp


template <class Tflow>
void CFlowNet<Tflow>::BI_Comp()
{
	int l,top,CompTop, CurNum=0;
	m_iCompNum=0;
	m_ipCompSep[0]=0;

	bool flag;

	for (int v=0; v < m_iVertNum; ++v) {
		m_ipCurEdge[v]=m_ipEdgeSep[v];
		m_ipNum[v]=NIL;
	}
	for (int e,v,w,i=0; i < m_iVertNum; ++i) {
		if (m_ipNum[i] == NIL) {
			if (m_ipEdgeSep[i+1]-m_ipEdgeSep[i] == 0) {
				m_ipComp[l=m_ipCompSep[m_iCompNum]]=i;
				m_ipCompSep[++m_iCompNum]=l+1;
				continue;
			}
			v=m_ipCompStack[CompTop=0]=m_ipStack[top=0]=i;
			m_ipBack[v]=m_ipNum[v]=CurNum++;
			for (; top >= 0; ) {
				l=m_ipEdgeSep[v+1];
				for (flag=true; m_ipCurEdge[v] < l;) {
					e=m_ipEdge[m_ipCurEdge[v]++];
					w=(m_ipTail[e] == v)? m_ipHead[e]: m_ipTail[e];						
					if (m_ipNum[w] == NIL) {
						v=m_ipCompStack[++CompTop]=m_ipStack[++top]=w;
						m_ipBack[v]=m_ipNum[v]=CurNum++;
						flag=false;
						break;
					}
					else {
						e=(top > 0)? m_ipStack[top-1]: NIL;
						if (e != w  &&  m_ipBack[v] > m_ipNum[w])
							m_ipBack[v]=m_ipNum[w];
					}
				}	// for (flag=true; m_ipCurEdge[v] < l;)
				if (flag && --top >= 0) {
					w=v;
					if (m_ipBack[v=m_ipStack[top]] > m_ipBack[w])
						m_ipBack[v]=m_ipBack[w];
					e=NIL;
					if (m_ipNum[w] == m_ipBack[w]) {
						flag=true;
						e=w;
					}
					else if (m_ipBack[w] == m_ipNum[v]) {
						flag=false;
						e=v;
					}
					if (e != NIL) {
						l=m_ipCompSep[m_iCompNum+1]=m_ipCompSep[m_iCompNum];
						for (w=m_ipCompStack[CompTop]; w != e; w=m_ipCompStack[CompTop]) {
							--CompTop;
							m_ipComp[l++]=w;
						}
						m_ipComp[l++]=w;
						if (flag)
							--CompTop;
						m_ipCompSep[++m_iCompNum]=l;
					}	// if (e != NIL)
				}	// if (flag && --top >= 0)
			}	// for (; top >= 0; )
		}	// if (Num[i] == NIL)
	}
}	// end of CFlowNet<Tflow>::BI_Comp

template <class Tflow>
void CFlowNet<Tflow>::BI_printSolution()
{
	int *Comp;
	int Size;
	cout << "Number of biconnected Components: " << m_iCompNum << endl;
	for (int i=0; i < m_iCompNum; ++i) {
		cout << i << ": ";
		Size=BI_getComp(i,Comp)-1;
		for (int j=0; j < Size; ++j)
			cout << Comp[j] << ", ";
		cout << Comp[Size] << endl;
	}	
}	// end of CFlowNet<Tflow>::BI_printSolution


//////////////////////////////////////////////////////////////
//             Minimal Cuts in undirected Graphs            //
//////////////////////////////////////////////////////////////
template <class Tflow>
void CFlowNet<Tflow>::printCut()
{

	int* parent=m_ipCurEdge;
	int i=0;
	cout << "Minimal Cut:\n";
	for (int v=0; v < m_iVertNum; ++v)
		if (parent[v] != NIL) {
			++i;
			cout << v;
			if ((i+1) % 5) cout << ", ";
			else cout << endl;
		}
	if ((i+1) % 5) cout << endl;
}	// CFlowNet<Tflow>::printCut

template <class Tflow>
void CFlowNet<Tflow>::MC_getMem()
{
	int iSize=sizeof(int)*(2*m_iMaxEdgeNum+4*m_iMaxVertNum+1)+
	(sizeof(Tflow)+sizeof(bool))*m_iMaxVertNum;
	if (!(m_pMC_MemBuf = new char[iSize])) {
		throw new CMemoryException("CFlowNet<Tflow>::MC_getMem");
	}
	m_ipEdgeSep=(m_ipEdge=reinterpret_cast<int*>(m_pMC_MemBuf))+(2*m_iMaxEdgeNum);
	m_ipNext=(m_ipName=m_ipEdgeSep+(m_iMaxVertNum+1))+m_iMaxVertNum;
	m_ipCut=m_ipNext+m_iMaxVertNum;
	m_ipCutCap=reinterpret_cast<Tflow*>(m_ipCut+m_iMaxVertNum);
	m_bpFlag=reinterpret_cast<bool*>(m_ipCutCap+m_iMaxVertNum);
}  // end of CFlowNet<Tflow>::MC_getMem

template <class Tflow>
void CFlowNet<Tflow>::MC_freeMem()
{
	if (m_pMC_MemBuf)
		delete[] (Tflow*)m_pMC_MemBuf;
}  // end of  CFlowNet<Tflow>::MC_freeMem

template <class Tflow>
void CFlowNet<Tflow>::MC_MinCut(Tflow threshold)
{
	for (int v=0; v < m_iVertNum; ++v) {
		m_ipName[v]=v;
		m_ipNext[v]=NIL;
	}

	int startVert=0, nextVert, e,k;
	Tflow maxCost;
	m_iCutVal=std::numeric_limits<Tflow>::max();

	for (int vertNum=m_iVertNum; vertNum > 1; --vertNum) { // main loop
		for (int v=0; v < m_iVertNum; ++v) {
			m_ipCutCap[e=m_ipName[v]]=0;
			m_bpFlag[e]=true;
		}
		m_bpFlag[startVert]=false;
		k=m_ipEdgeSep[startVert+1];
		for (int v=m_ipEdgeSep[startVert]; v < k; ++v) {
			e=m_ipEdge[v];
			m_ipCutCap[m_ipName[getOtherEnd(e,startVert)]]+=m_ipUCap[e];
		}

		for (int i=2; i < vertNum; ++i) {
			// looking for the next vertex
			maxCost=-std::numeric_limits<Tflow>::max();

			for (int v=0; v < m_iVertNum; ++v) {
				e=m_ipName[v];
				if (m_bpFlag[e] && (m_ipCutCap[e] > maxCost)) {
					maxCost = m_ipCutCap[e];
					nextVert=e;
				}
			}
			m_bpFlag[nextVert]=false;
			// update CuTflow 
			for (int v=nextVert; v >= 0; v = m_ipNext[v]) {
				k=m_ipEdgeSep[v+1];
				for (int i=m_ipEdgeSep[v]; i < k; ++i) {
					e=m_ipEdge[i];
					m_ipCutCap[m_ipName[getOtherEnd(e,v)]]+=m_ipUCap[e];
				}
			}
		}	// for (i=2; i < vertNum; ++i)

		for (int v=0; v < m_iVertNum; ++v)
			if (m_bpFlag[m_ipName[v]]) {
				k=v; break;
			}
		if (m_ipCutCap[k=m_ipName[k]] < m_iCutVal) {
			m_iCutVal=m_ipCutCap[k];
			memset(m_ipCut,0,m_iVertNum*sizeof(int));
			for (int i=k; i >= 0; i = m_ipNext[i])
				m_ipCut[i]=1;
			if (isPositive(threshold-m_iCutVal))
				return;
		}
		// shrink two last vertices in the sequence
		for (int i=k; i >=0; i = m_ipNext[i]) {
			m_ipName[i]=nextVert;
		}
		for (int i=nextVert; i >= 0; i = m_ipNext[e=i]);
		m_ipNext[e]=k;
	}  // end of the main loop
}  // end of CFlowNet<Tflow>::MC_MinCut

template <class Tflow>
void CFlowNet<Tflow>::MC_printSolution()
{
	cout << "\nMinimum Cut Value = " << m_iCutVal << endl;
	cout << "\nMinimum Cut:\n(";
	int n=m_iVertNum-1;
	for (int i=0; i < n; i++)
		cout << m_ipCut[i] << ",";
	cout << m_ipCut[n] << ")\n";
}	// end of CFlowNet<Tflow>::MC_printSolution

/////////////////////////////////////////////////////////
//                     F L O W S                       //
/////////////////////////////////////////////////////////
template <class Tflow>
Tflow CFlowNet<Tflow>::getLoCap(int e)
{
	return
	(m_iProblemMsk & m_iDirectMsk)?
			((m_iProblemMsk & m_iLoCapMsk)? m_ipLCap[e]: 0): -m_ipUCap[e];
}	// end of CFlowNet<Tflow>::getLoCap


template <class Tflow>
Tflow CFlowNet<Tflow>::resCap(int e, int v)
{
	return
	(v==m_ipTail[e])? m_ipUCap[e]-m_ipFlow[e]: m_ipFlow[e]-getLoCap(e);
}  // end of CFlowNet<Tflow>::resCap

template <class Tflow>
void CFlowNet<Tflow>::printFlow()
{
	cout << "Flow:\n";
	int i;
	for (i=0; i < m_iEdgeNum; i++) {
		cout << "f(" << m_ipTail[i] << "," << m_ipHead[i]
		                                               << ")=" << m_ipFlow[i];
		if ((i+1) % 5) cout << ", ";
		else cout << endl;
	}
	if ((i+1) % 5) cout << endl;
}	// end of CFlowNet<Tflow>::PrintFlow


/////////////////////////////////////////////////////////
//                  PUSH and RELABEL                   //
/////////////////////////////////////////////////////////

template <class Tflow>
void CFlowNet<Tflow>::PR_push(int v, int e)
// returns true in case of saturated Push
{
	Tflow delta=min(resCap(e,v),m_ipExcess[v]);
	int w=getOtherEnd(e,v);

	if (v==m_ipTail[e])
		m_ipFlow[e]+=delta;
	else
		m_ipFlow[e]-=delta;

	m_ipExcess[v]-=delta;

	if (w != m_iSource && w!= m_iThink && isZero(m_ipExcess[w]))
		addToQueue(w);
	m_ipExcess[w]+=delta;
}  // end of CFlowNet<Tflow>::PR_push

template <class Tflow>
void CFlowNet<Tflow>::PR_relabel(int v)
{
	int wl=std::numeric_limits<int>::max();
	int e, k=m_ipEdgeSep[v+1];
	for (int i=m_ipEdgeSep[v]; i < k; i++)
		if (isPositive(resCap(e=m_ipEdge[i],v)))
			wl=min(wl,m_ipLabel[getOtherEnd(e,v)]);

	m_ipLabel[v]=wl+1;
}  // end of CFlowNet<Tflow>::PR_relabel

template <class Tflow>
int CFlowNet<Tflow>::PR_getNextFeasible(int v, int &i)
{
	int p=m_ipEdgeSep[v];
	int q=m_ipEdgeSep[v+1];
	int w,e,e1=NIL;
	for (int k=p; k < q; ++k) {
		w=getOtherEnd(e=m_ipEdge[i],v);
		if (isPositive(resCap(e,v)) && (m_ipLabel[v] == m_ipLabel[w] + 1)) {
			e1=e;
			break;
		}
		if (++i == q) i=p;
	}
	return e1;
}  // end of CFlowNet<Tflow>::PR_getNextFeasible

template <class Tflow>
void CFlowNet<Tflow>::PR_discharge(int v)
{
	int p=m_ipEdgeSep[v+1];
	int q=m_ipEdgeSep[v];
	int i;
	bool flag=true;
	int NoOfOutArcsInGf=0;
	int e;
	for (i=m_ipCurEdge[v]; ;) {
		if ((e=PR_getNextFeasible(v,i)) == NIL) {
			PR_relabel(v);
			addToQueue(v);
			break;
		}
		PR_push(v,e);
		if (isZero(m_ipExcess[v]))
			break;
	}
	m_ipCurEdge[v]=i; 
}  // end of CFlowNet<Tflow>::PR_discharge

template <class Tflow>
bool CFlowNet<Tflow>::findFeasibleSol()
{
	if (!(m_iProblemMsk & m_iLoCapMsk)) {
		memset(m_ipFlow,0,m_iEdgeNum*sizeof(Tflow));
		return true;
	}
	int q;
	Tflow* upCap=m_ipUCap+m_iEdgeNum;
	memcpy(upCap,m_ipDemand,m_iVertNum*sizeof(Tflow));
	memset(m_ipLCap+m_iEdgeNum,0,m_iVertNum*sizeof(Tflow));

	for (int e=0; e < m_iEdgeNum; ++e) {
		upCap[m_ipHead[e]]-=(m_ipFlow[e]=m_ipLCap[e]);
		upCap[m_ipTail[e]]+=m_ipFlow[e];
	}

	m_iSource=m_iVertNum;
	m_iThink=m_iVertNum+1;

	for (int e=0; e < m_iVertNum; ++e) {
		m_ipFlow[q=m_iEdgeNum+e]=0;
		if (isNonNegative(upCap[e])) {
			m_ipTail[q]=e;
			m_ipHead[q]=m_iThink;
		}
		else {
			upCap[e]=-upCap[e];
			m_ipTail[q]=m_iSource;
			m_ipHead[q]=e;
		}
	}
	m_iEdgeNum+=m_iVertNum;
	m_iVertNum+=2;
	buildEdgeList();
	PR_pushAndRelabel();
	q=PR_BDF();

	m_iVertNum-=2;
	m_iEdgeNum-=m_iVertNum;
	return (q-1)? false: true;
}  // end of CFlowNet<Tflow>::FindFeasibleSol

template <class Tflow>
void CFlowNet<Tflow>::PR_initSolution()
{
	memset(m_ipLabel,0,m_iVertNum*sizeof(int));
	m_ipLabel[m_iSource]=m_iVertNum;
	memset(m_ipExcess,0,m_iVertNum*sizeof(Tflow));

	memcpy(m_ipCurEdge,m_ipEdgeSep,m_iVertNum*sizeof(int));

	initQueue();

	int e,v, k=m_ipEdgeSep[m_iSource+1];
	Tflow uf,LoCap;
	for (int i=m_ipEdgeSep[m_iSource]; i < k; i++)
		if (m_ipTail[e=m_ipEdge[i]] == m_iSource) {
			if (isPositive(uf=m_ipUCap[e]-m_ipFlow[e])) {
				m_ipExcess[v=m_ipHead[e]]+=uf;
				m_ipExcess[m_iSource]-=uf;
				m_ipFlow[e]=m_ipUCap[e];
				if (v !=m_iThink)
					addToQueue(v);
			}
		}
		else { // m_ipHead[e] != m_iSource
			if (isPositive(uf=m_ipFlow[e]-(LoCap=getLoCap(e)))) {
				m_ipExcess[v=m_ipTail[e]]+=uf;
				m_ipExcess[m_iSource]-=uf;
				m_ipFlow[e]=LoCap;
				if (v !=m_iThink)
					addToQueue(v);
			}
		}
} // end of CFlowNet<Tflow>::PR_initSolution

template <class Tflow>
void CFlowNet<Tflow>::initQueue()
{
	m_iFirst=m_iLast=0;
}

template <class Tflow>
void CFlowNet<Tflow>::addToQueue(int v)
{
	m_ipQueue[m_iLast]=v;
	m_iLast=(m_iLast+1) % m_iVertNum;
}  // end of CFlowNet<Tflow>::AddToQueue

template <class Tflow>
int CFlowNet<Tflow>::getFromQueue()
{
	int v=(m_iFirst==m_iLast)? NIL: m_ipQueue[m_iFirst];
	m_iFirst=(m_iFirst+1) % m_iVertNum;
	return v;  
}  // end of CFlowNet<Tflow>::GetFromQueue

template <class Tflow>
void CFlowNet<Tflow>::PR_getMem()
{
	int iSize=sizeof(int)*(2*m_iMaxEdgeNum+4*m_iMaxVertNum+1)+
	sizeof(Tflow)*(m_iMaxEdgeNum+m_iMaxVertNum);
	if (!(m_pPR_MemBuf = new int[iSize])) {
		throw new CMemoryException("CFlowNet<Tflow>::PR_getMem");
	}

	m_ipEdgeSep=(m_ipEdge=(int*)m_pPR_MemBuf)+2*m_iMaxEdgeNum;
	m_ipQueue=(m_ipCurEdge=(m_ipLabel=m_ipEdgeSep+(m_iMaxVertNum+1))+m_iMaxVertNum)+m_iMaxVertNum;
	m_ipExcess=(m_ipFlow=(Tflow*)(m_ipQueue+m_iMaxVertNum))+m_iMaxEdgeNum;
}  // end of  CFlowNet<Tflow>::PR_getMem

template <class Tflow>
void CFlowNet<Tflow>::PR_freeMem()
{
	if (m_pPR_MemBuf)
		delete[] (Tflow*)m_pPR_MemBuf;
}  // end of  CFlowNet<Tflow>::PR_freeMem

template <class Tflow>
int CFlowNet<Tflow>::PR_BDF()
{
	int* Parent=m_ipCurEdge;
	int e,w;
	initQueue();
	addToQueue(m_iSource);
	for (int v=0; v < m_iVertNum; ++v)
		Parent[v]=NIL;
	Parent[m_iSource]=m_iSource;
	int Size=1;
	for (int v=getFromQueue(); v != NIL; v=getFromQueue())
		for (int i=m_ipEdgeSep[v]; i < m_ipEdgeSep[v+1]; i++)
			if (isPositive(resCap(e=m_ipEdge[i],v)) && Parent[w=getOtherEnd(e,v)] == NIL) {
				//      if (resCap(e=m_ipEdge[i],v) > 0 && Parent[w=GetOtherEnd(e,v)] == NIL) {
				addToQueue(w);
				Parent[w]=v;
				Size++;
			}
	return Size;
}  // end of CFlowNet<Tflow>::PR_BDF

template <class Tflow>
void CFlowNet<Tflow>::PR_pushAndRelabel()
{
	PR_initSolution();
	for (int v=getFromQueue(); v != NIL; v=getFromQueue())
		PR_discharge(v);
	m_iMaxFlow=m_ipExcess[m_iThink];
}  // end of CFlowNet<Tflow>::PushAndRelabel

template <class Tflow>
void CFlowNet<Tflow>::PR_printSolution()
{
	printFlow();
	printCut();
}  // end of CFlowNet<Tflow>::PR_printSolution


///////////////////////////////////////////////////
//        Homory Cut Tree: Implementation        //
///////////////////////////////////////////////////
template <class Tflow>
void CFlowNet<Tflow>::GH_getMem()
{
	int iSize=(sizeof(int)+sizeof(Tflow))*m_iMaxVertNum;
	if (!(m_pGH_MemBuf = new char[iSize]))
		throw new CMemoryException(
				"CFlowNet<Tflow>::MC_getMem: Error allocating memory");

	m_ipFl=(Tflow*)((m_ipP = (int*)m_pGH_MemBuf)+m_iMaxVertNum);
}	// end of CFlowNet<Tflow>::GH_getMem

template <class Tflow>
void CFlowNet<Tflow>::GH_freeMem()
{
	if (m_pGH_MemBuf)
		delete[] (Tflow*)m_pGH_MemBuf;
}	// end of CFlowNet<Tflow>::GH_freeMem

template <class Tflow>
void CFlowNet<Tflow>::GH_GomoryTree()
{
	int* Parent=m_ipCurEdge;
	int t;
	memset(m_ipP,0,sizeof(int)*m_iVertNum);
	for (int s=1; s < m_iVertNum; ++s) {
		memset(m_ipFlow,0,sizeof(Tflow)*m_iEdgeNum);
		t=m_iThink=m_ipP[m_iSource=s];
		PR_pushAndRelabel();
		(void)PR_BDF();
		for (int i=0; i < m_iVertNum; ++i)
			if (i != s && Parent[i] != NIL && m_ipP[i] == t)
				m_ipP[i]=s;

		if (Parent[m_ipP[t]] != NIL) {
			m_ipP[s]=m_ipP[t];
			m_ipP[t]=s;
			m_ipFl[s]=m_ipFl[t];
			m_ipFl[t]=m_iMaxFlow;
		}	// if (Parent[m_ipP[t]] != NIL)
		else
			m_ipFl[s]=m_iMaxFlow;

	}	// for (register int s=1;	
	//_ipP[0]=NIL;	// to mark the root
}	// end of CFlowNet<Tflow>::GH_GomoryTree


template <class Tflow>
void CFlowNet<Tflow>::GH_printSolution()
{
	cout << "Homory Cut Tree:\n";
	for (int i=1; i < m_iVertNum; i++) {
		cout << "fl(" << i << "," << m_ipP[i] << ")=" << m_ipFl[i];
		if (!((i+1) % 5))
			cout << endl;
		else
			cout << ", ";
	}
	cout << endl;
}	// end of GH_printSolution

///////////////////////////////////////////////////
//          Odd Cut: Implementation              //
///////////////////////////////////////////////////

template <class Tflow>
inline void CFlowNet<Tflow>::DS_MakeSet(int x)
{
	m_ipPrec[x]=x; m_ipDepth[x]=0;
}	// end of CFlowNet<Tflow>::DS_MakeSet

template <class Tflow>
int CFlowNet<Tflow>::DS_Find(int x)
// this program implement heuristic, called path halving
{
	while (m_ipPrec[m_ipPrec[x]] != m_ipPrec[x])
		x=m_ipPrec[x]=m_ipPrec[m_ipPrec[x]];
	return m_ipPrec[x];

}	// end of CFlowNet<Tflow>::DS_Find

template <class Tflow>
int CFlowNet<Tflow>::DS_Link(int x, int y)
// this program implement heuristic, called union by rank
{
	if (m_ipDepth[x] > m_ipDepth[y]) {
		int z=x;
		x=y; y=z;
	}
	else if (m_ipDepth[x] == m_ipDepth[y]) m_ipDepth[y]++;
	return (m_ipPrec[x]=y);

}	// end of CFlowNet<Tflow>::DS_Link


template <class Tflow>
void CFlowNet<Tflow>::TCUT_getMem()
{
	try {
		PR_getMem();
		GH_getMem();
	}
	catch(CMemoryException* pe) {
		throw pe;
	}
}	// end of TCUT_getMem

template <class Tflow>
void CFlowNet<Tflow>::TCUT_freeMem()
{
	PR_freeMem();
	GH_freeMem();
}	// end of TCUT_freeMem

template <class Tflow>
void CFlowNet<Tflow>::T_Cut()
{
	GH_GomoryTree();
	//	GH_printSolution();
	m_ipPrec=m_ipCurEdge;
	m_ipDepth=m_ipLabel;
	int Card,j,side;
	m_ipCut=reinterpret_cast<int*>(m_ipExcess);
	m_iCutVal=std::numeric_limits<Tflow>::max();
	int ir;
	for (int i=1; i < m_iVertNum; ++i)
		if (isPositive(m_iCutVal - m_ipFl[i])) {
			for (int v=0; v < m_iVertNum; ++v) {
				DS_MakeSet(v);
			}
			for (int v=1; v < m_iVertNum; ++v) {
				if (v != i)
					DS_Link(DS_Find(v),DS_Find(m_ipP[v]));
			}

			ir=DS_Find(i);
			for (int v=Card=0; v < m_iTermNum; ++v)
				if (ir == DS_Find(m_ipTerminal[v]))
					Card++;
			if (Card % 2) {
				side=(2*Card < m_iVertNum)? 1: 0;
				m_iCutVal=m_ipFl[j=i];
			}
		}	// if (BestCut < m_ipFl[i])

	for (int i=0; i < m_iVertNum; ++i) {
		DS_MakeSet(i);
	}
	for (int i=1; i < m_iVertNum; ++i) {
		if (i != j)
			DS_Link(DS_Find(i),DS_Find(m_ipP[i]));
	}

	for (int i=0; i < m_iVertNum; ++i) {
		m_ipCut[i]=(DS_Find(j) == DS_Find(i))? side: 1-side;
	}
}	// end of CFlowNet<Tflow>::T_Cut()


template <class Tflow>
void CFlowNet<Tflow>::COMB_getMem()
{
	try {
		PR_getMem();
		GH_getMem();
	}
	catch(CMemoryException* pe) {
		throw pe;
	}
}	// end of COMB_getMem

template <class Tflow>
void CFlowNet<Tflow>::COMB_freeMem()
{
	PR_freeMem();
	GH_freeMem();
}	// end of COMB_freeMem


template <class Tflow>
void CFlowNet<Tflow>::COMB_AND_CUT_getMem()
{
	int iSize=sizeof(int)*(2*m_iMaxEdgeNum+6*m_iMaxVertNum+1)+
	(2*m_iMaxVertNum+m_iMaxEdgeNum)*sizeof(Tflow)+ 
	m_iMaxVertNum*(sizeof(bool));
	if (!(m_pMC_MemBuf = new char[iSize])) {
		throw new CMemoryException("CFlowNet<Tflow>::COMB_AND_CUT_getMem");
	}
	m_ipEdgeSep=(m_ipEdge=(int*)m_pMC_MemBuf)+(2*m_iMaxEdgeNum);
	m_ipNext=(m_ipName=m_ipEdgeSep+(m_iMaxVertNum+1))+m_iMaxVertNum;
	m_ipP=(m_ipQueue=m_ipNext+m_iMaxVertNum)+m_iMaxVertNum;
	m_ipCut=m_ipP+m_iMaxVertNum;

	m_ipFl=(m_ipFlow=(m_ipCutCap=(Tflow*)(m_ipCut+m_iMaxVertNum))+m_iMaxVertNum)+m_iMaxEdgeNum;

	m_ipLabel=m_ipName;
	m_ipCurEdge=m_ipNext;
	m_ipExcess=m_ipCutCap;

	m_bpFlag=(bool*)(m_ipFl+m_iMaxVertNum);
}	// end of COMB_AND_CUT_getMem()

////////////////////////////////////////////////////////////////////
// G^*(x)=(V^*(x),E^*(x)), V^*(x) = V \cup {v_e: e\in E},         //
//        E^*(x) = {(v,v_e),(v_e,w): e=(v,w) in E, x(e) > 0}      //
//                     u(v,v_e) = x(e), u(v_e,w) = 1-x(e)         //
// Odd = {v_e: e in E} \cup {w in V: |{(v_e,w): e in E}| is odd}  //
//----------------------------------------------------------------//
// Algorithm:                                                     //
// 1) find a minimum odd cut (S,V\S) in (G^*(x),u)                //
//    (i.e., such that |S\cap Odd| is odd)                        //
// 2) build the comb (H; T_1,...,T_{2K+1}:                        //
//    a) all nodes from S \cap V belong to the handle H;          //
//    b) at most one of the two edges (v,v_e) and (v_e,w) belongs //
//       to the cut (S,V\S), if it is (v_e,w), then e=(v,w)       //
//       is taken as a tooth;                                     //
//    c) if two teeth intersect in a common node q,               //
//       the both teeth are removed from the structure;           //
//       in addition, if q is in the handle H,                    //
//       then it is removed from it; otherwise, it is added to it //
////////////////////////////////////////////////////////////////////

template <class Tflow>
int CFlowNet<Tflow>::getNextBlossom(int &s, double delta,
		int &hdSize, int* &ipComb)
{
	int iEdgeNum=m_iEdgeNum/2;
	int	iVertNum=m_iVertNum-iEdgeNum;
	m_ipPrec=m_ipCurEdge;
	m_ipDepth=m_ipLabel;
	int e,ir;
	ipComb=(int*)m_ipExcess;
	bool* bpTerminal=(bool*)m_ipFlow;

	int i,v;
	for (i=s; i < m_iVertNum; i++) {
		if (1- m_ipFl[i] > delta) {
			for (v=0; v < m_iVertNum; v++)
				DS_MakeSet(v);
			for (v=1; v < m_iVertNum; v++)
				if (v != i) DS_Link(DS_Find(v),DS_Find(m_ipP[v]));

			ir=DS_Find(i);

			for (e=v=0; v < iVertNum; v++) {
				if (bpTerminal[v])
					if (ir == DS_Find(v))
						e++;
			}

			for (v=iVertNum; v < m_iVertNum; v++)
				if (ir == DS_Find(v))
					++e;

			if (e % 2)
				break;
		}
	}
	s=i+1;
	if (i >= m_iVertNum) {
		for (i=0,e=iEdgeNum; i < iEdgeNum;)
			m_ipHead[i++]=m_ipHead[e++];
		m_iVertNum=iVertNum; m_iEdgeNum=iEdgeNum;
		return 0;
	}

	int w, z=0;

	for (i=0; i < iVertNum; ++i)
		ipComb[i]=(ir == DS_Find(i))? m_iVertNum: -m_iVertNum;

	for (e=0, i=iEdgeNum; i < m_iEdgeNum; ++e,++i) {
		if (DS_Find(m_ipTail[i]) != DS_Find(m_ipHead[i])) {
			if ((z=abs(ipComb[v=m_ipTail[e]])-1) < iVertNum || 
					abs(ipComb[w=m_ipHead[i]]) <= iVertNum) {
				if (z >= iVertNum)
					z=abs(ipComb[v=w])-1;
				ipComb[z]=ipComb[v]=(ipComb[v] < 0)? m_iVertNum: -m_iVertNum;
			}
			else {
				ipComb[v]=(ipComb[v] > 0)? w+1: -w-1;
				ipComb[w]=(ipComb[w] > 0)? v+1: -v-1;
			}
		}
	}	// for (e=0;...

	for (w=i=0; i < iVertNum; ++i) {
		if (ipComb[i] > 0) {
			ipComb[i]=1;
			w++;
		}
	}
	hdSize=w;
	for (z=i=0; i < iVertNum; ++i) // z - counter for teeth
		if ((v=-ipComb[i]) > 0) {
			if ((--v) < iVertNum) {
				++z;
				ipComb[i]=z << 16;
				ipComb[v]|=(z << 16);
			}
			else ipComb[i]=0;
		}

	return z;
}	// end of getNextBlossom

template <class Tflow>
void CFlowNet<Tflow>::blossom()
{
	int v=m_iVertNum, e=m_iEdgeNum;
	for (int i=0; i < m_iEdgeNum;) {
		m_ipUCap[e]=1.0-m_ipUCap[i];
		m_ipHead[e]=m_ipHead[i];
		m_ipHead[i++]=m_ipTail[e++]=v++;
	}
	int iVertNum=m_iVertNum, iEdgeNum=m_iEdgeNum;
	m_iVertNum=v; m_iEdgeNum=e;
	buildEdgeList();
	GH_GomoryTree();
	//	GH_printSolution();

	bool* bpTerminal=(bool*)m_ipFlow;

	for (e=0; e < iVertNum; e++)
		bpTerminal[e]=false;

	for (e=iEdgeNum; e < m_iEdgeNum; e++) {
		v=m_ipHead[e];
		bpTerminal[v]=!bpTerminal[v];
	}

}	// end of blossom

// computes maximally violated (i.e. by 1) blossom inequalities
template <class Tflow>
int CFlowNet<Tflow>::maxBlossom()
{
	m_ipPrec=m_ipCurEdge;
	m_ipDepth=m_ipLabel;
	int *ipNext=m_ipP, *ipEdge=(int*)m_ipFl, *ipDeg=(int*)m_ipFlow;

	for (int i=0; i < m_iVertNum; ++i) {
		ipNext[i]=i;
		DS_MakeSet(i);
	}

	int e,k,r1, r2, v,w;
	memset(ipDeg,0,m_iVertNum*sizeof(int));
	for (int i=k=0; i < m_iEdgeNum; ++i) {
		if (isZero(1-m_ipUCap[i])) {
			ipDeg[m_ipTail[i]]++;
			ipDeg[m_ipHead[i]]++;
			ipEdge[k++]=i;
		}
		else if ((r1=DS_Find(m_ipTail[i])) != (r2=DS_Find(m_ipHead[i])))
			DS_Link(r1,r2);
	}
	for (int i=v=0; i < k; ++i) {
		e=ipEdge[i];
		if (ipDeg[m_ipTail[e]] > 1 && ipDeg[m_ipHead[e]] > 1) {
			if ((r1=DS_Find(m_ipTail[e])) != (r2=DS_Find(m_ipHead[e])))
				DS_Link(r1,r2);
		}
		else ipEdge[v++]=e;
	}
	k=v;
	for (int i=0; i < k; ++i) {
		e=ipEdge[i];
		r1=DS_Find(m_ipTail[e]);
		r2=DS_Find(m_ipHead[e]);
		if (ipDeg[r1] > ipDeg[r2]) {
			v=r1; r1=r2; r2=v;
		}
		if (ipDeg[r2] == 2)
			ipDeg[r2]=r1+3;
		else if (ipDeg[r2]=r1+3)
			DS_Link(r1,r2);
	}	
	memset(ipDeg,0,m_iVertNum*sizeof(int));

	for (int i=0; i < k; ++i) {
		e=ipEdge[i];
		++ipDeg[DS_Find(m_ipTail[e])];
		++ipDeg[DS_Find(m_ipHead[e])];
	}

	int* ipCombs=(int*)m_ipFl;
	for (int i=0; i < m_iVertNum; ++i) {
		ipCombs[i]=-1;
		if ((r1=DS_Find(i)) != i) {
			ipNext[i]=ipNext[r1];
			ipNext[r1]=i;
		}
	}

	for (int i=0; i < m_iEdgeNum; ++i) {
		if (isZero(1-m_ipUCap[i]))
			if ((r1=DS_Find(v=m_ipTail[i])) != (r2=DS_Find(w=m_ipHead[i]))) {
				if (ipDeg[r1] % 2)
					ipCombs[v]=w;
				if (ipDeg[r2] % 2)
					ipCombs[w]=v;
			}
	}
	for (int i=r1=0; i < m_iVertNum; ++i) {
		if (ipDeg[i] % 2)
			m_ipDepth[r1++]=i;
	}

	return r1;
}	// end of maxBlossom

template <class Tflow>
void CFlowNet<Tflow>::getMaxBlossom(int k, int &hdSize, int &thNum, int* &ipComb)
{
	int *ipNext=m_ipP, *ipCombs=(int*)m_ipFl, *ipRoot=m_ipLabel;
	hdSize=thNum=0;
	ipComb=(int*)m_ipExcess;
	memset(ipComb,0,m_iVertNum*sizeof(int));
	for (int i=k=ipRoot[k];;) {
		ipComb[i]=1; hdSize++;
		if (ipCombs[i] >= 0) {
			++thNum;
			ipComb[ipCombs[i]]=(thNum << 16);
			ipComb[i]|=(thNum << 16);
		}
		if ((i=ipNext[i]) == k)
			break;
	}
}	// end of getMaxBlossom(int k)

template <class Tflow>
void CFlowNet<Tflow>::COMB_printSolution()
{
	int* ipComb=(int*)m_ipExcess;
	int k;
	int p=m_iVertNum/2+1;
	cout << "Handle:\n";
	for (int i=k=0; i < m_iVertNum; i++)
		if (ipComb[i] >= p) {
			cout << i;
			if (!(k=(++k % 10)))
				cout << endl;
			else
				cout << ", ";
		}

	cout << "\nTeeth:\n";
	for (int j,i=k=0; i < m_iVertNum; i++)
		if (j=ipComb[i] % p) {
			cout << "(" << i << " in " << j << ")";
			if (!(k=(++k % 10)))
				cout << endl;
			else
				cout << ", ";
		}
	cout << endl;
}	// end of COMB_printSolution


/////////////////////////////////////////////
//         Network Simplex Method          //
/////////////////////////////////////////////

template <class Tflow>
void CFlowNet<Tflow>::TRAN_getMem()
{
	int iSize=sizeof(Tflow)*(m_iMaxEdgeNum+m_iMaxVertNum)+
	sizeof(int)*3*m_iMaxVertNum;
	if (!(m_pTRAN_MemBuf=new char[iSize]))
		throw new CMemoryException(
				"CFlowNet<Tflow>::TRAN_getMem: Error allocating memory");

	m_ipPrice=(m_ipFlow=(Tflow*)m_pTRAN_MemBuf)+m_iMaxEdgeNum;
	m_ipThread=(m_ipDepth=(m_ipPrec=(int*)(m_ipPrice+m_iMaxVertNum))+m_iMaxVertNum)+m_iMaxVertNum;
}	// end of CFlowNet<Tflow>::TRAN_getMem

template <class Tflow>
void CFlowNet<Tflow>::TRAN_freeMem()
{
	if (m_pTRAN_MemBuf)
		delete[] (Tflow*)m_pTRAN_MemBuf;
}	// end of CFlowNet<Tflow>::TRAN_freeMem

template <class Tflow>
int CFlowNet<Tflow>::TRAN_delEdge(int edge, int delta)
{

	int v, w, x, y;
	if (m_ipPrec[w=m_ipHead[edge]] == edge)
		v=m_ipTail[edge];
	else {
		v=w;
		w=m_ipTail[edge];
	}
	int gamma=m_ipDepth[w];
	m_ipDepth[w]=delta;
	delta-=gamma;
	for (x=w; m_ipDepth[y=m_ipThread[x]] > gamma; x=y)
		m_ipDepth[y]+=delta;
	m_ipThread[x]=w;
	m_ipPrec[w]=NIL;
	// looking for the last vertex in the subtree
	// rooted at v
	for (; m_ipThread[v] != w; v=m_ipThread[v]);
	m_ipThread[v]=y;
	return x;
}	// end of CFlowNet<Tflow>::TRAN_delEdge

template <class Tflow>
void CFlowNet<Tflow>::TRAN_addEdge(int edge, int direct)
{
	int v, w, x, y, z;
	if (direct == 1) {
		v=m_ipTail[edge];
		w=m_ipHead[edge];
	}
	else {
		v=m_ipHead[edge];
		w=m_ipTail[edge];
	}
	int h=m_ipDepth[v];

	for (x=v; m_ipDepth[y=m_ipThread[x]] > h; x=y);

	int edge1;
	for (;(edge1=m_ipPrec[w]) != NIL;) {
		z=TRAN_delEdge(edge1,m_ipDepth[v]+1);
		m_ipPrec[w]=edge;
		m_ipThread[x]=v=w;
		x=z;
		edge=edge1;
		w=(w==m_ipTail[edge])? m_ipHead[edge]: m_ipTail[edge];
	}
	m_ipPrec[w]=edge;
	m_ipThread[x]=w;
	h=m_ipDepth[v]+1;
	for (x=w; m_ipThread[x] != w; x=m_ipThread[x])
		m_ipDepth[x]+=h;
	m_ipThread[x]=y;
	m_ipDepth[x]+=h;
}	// end of CFlowNet<Tflow>::TRAN_addEdge

template <class Tflow>
int CFlowNet<Tflow>::TRAN_findRoot(int vertex)
{
	int e;
	for (; m_ipPrec[vertex] != NIL;) {
		e=m_ipPrec[vertex];
		vertex=(m_ipTail[e])? m_ipHead[e]: m_ipTail[e];
	}
	return vertex;
}	// end of CFlowNet<Tflow>::TRAN_findRoot

template <class Tflow>
void CFlowNet<Tflow>::getEmptyTree()
{
	for (int v=0; v < m_iVertNum; v++) {
		m_ipPrec[v]=NIL;
		m_ipDepth[v]=0;
		m_ipThread[v]=v;
	}
}	// end of CFlowNet<Tflow>::getEmptyTree

template <class Tflow>
void CFlowNet<Tflow>::TRAN_initTree()
{
	getEmptyTree();
	for (int e=0; e < m_iEdgeNum; ++e) {
		if (TRAN_findRoot(m_ipTail[e]) != TRAN_findRoot(m_ipHead[e]))
			TRAN_addEdge(e,1);
	}
}	// end of CFlowNet<Tflow>::TRAN_initTree

template <class Tflow>
void CFlowNet<Tflow>::TRAN_computePrices()
{
	for (int e,v=0; v < m_iVertNum; ++v) {
		if (m_ipPrec[v] == NIL) {
			m_ipPrice[v]=0;
			for (int w=m_ipThread[v]; w != v; w=m_ipThread[w]) {
				e=m_ipPrec[w];
				m_ipPrice[w] = (w == m_ipHead[e])?
						m_ipPrice[m_ipTail[e]] + m_ipCost[e]:
				m_ipPrice[m_ipHead[e]] - m_ipCost[e];
			}
		}
	}
}	// end of CFlowNet<Tflow>::TRAN_computePrices

template <class Tflow>
int CFlowNet<Tflow>::TRAN_changeFlow(int edge, int& side, int direct)
{
	int v,w,x,y; 
	int e1, e2, edg;
	Tflow delta,delta1,delta2;
	side=1;
	if (direct == 1) {
		v=m_ipTail[edge];
		w=m_ipHead[edge];
	}
	else {
		v=m_ipHead[edge];
		w=m_ipTail[edge];
	}
	delta1=resCap(edge,v);

	delta2=delta1; e1=edge; e2=edge;
	x=v; y=w;

	for (; x!=y;)
		if (m_ipDepth[x] >= m_ipDepth[y]) {
			if (x == m_ipHead[edg=m_ipPrec[x]]) {
				delta=m_ipUCap[edg]-m_ipFlow[edg];
				x=m_ipTail[edg];
			}
			else {
				delta=m_ipFlow[edg]-getLoCap(edg);
				x=m_ipHead[edg];
			}
			if (isNegative(delta-delta1)) {
				delta1=delta;
				e1=edg;
			}
		}
		else {
			if (y == m_ipTail[edg=m_ipPrec[y]]) {
				delta=m_ipUCap[edg]-m_ipFlow[edg];	
				y=m_ipHead[edg];	
			}
			else {
				delta=m_ipFlow[edg]-getLoCap(edg);
				y=m_ipTail[edg];
			}
			if (isNegative(delta-delta2)) {
				delta2=delta;
				e2=edg;
			}
		}
	if (isPositive(delta1-delta2)) {
		delta1=delta2;
		e1=e2;
		side=-1;
	}
	m_ipFlow[edge]+=(direct == 1)? delta1: -delta1;
	for (x=v; x != y; ) {
		edg=m_ipPrec[x];
		if (x == m_ipHead[edg]) {
			m_ipFlow[edg]+=delta1;
			x=m_ipTail[edg];
		}
		else {
			m_ipFlow[edg]-=delta1;
			x=m_ipHead[edg];
		}
	}
	for (x=w; x != y; ) {
		edg=m_ipPrec[x];
		if (x == m_ipTail[edg]) {
			m_ipFlow[edg]+=delta1;
			x=m_ipHead[edg];
		}
		else {
			m_ipFlow[edg]-=delta1;
			x=m_ipTail[edg];
		}
	}
	return e1;
}	// end of CFlowNet<Tflow>::TRAN_changeFlow

template <class Tflow>
int CFlowNet<Tflow>::TRAN_checkFlowOpt(int& direct)
{
	int v,w;
	Tflow cp;
	for (int e=0; e < m_iEdgeNum; ++e) {
		if (e == m_ipPrec[m_ipTail[e]] || e == m_ipPrec[m_ipHead[e]])
			continue;
		cp=m_ipPrice[v=m_ipTail[e]] + m_ipCost[e] - m_ipPrice[w=m_ipHead[e]];
		if (isNegative(cp) && isPositive(m_ipUCap[e]-m_ipFlow[e])) {
			direct=1;
			return e;
		}
		if (isPositive(cp) && isPositive(m_ipFlow[e]-m_ipLCap[e])) {
			direct=-1;
			return e;
		}
	}
	return NIL;
}	// CFlowNet<Tflow>::TRAN_checkFlowOpt

template <class Tflow>
void CFlowNet<Tflow>::TRAN_setBigM()
{
	Tflow delta=1;
	for (int e=0; e < m_iEdgeNum; ++e)
		delta+=abs(m_ipCost[e]);
	m_iBigM=delta;
}	// end of CFlowNet<Tflow>::TRAN_setBigM

template <class Tflow>
void CFlowNet<Tflow>::TRAN_initSolution()
{
	Tflow* Excess=(Tflow*)m_ipPrice;
	Tflow delta;

	memcpy(Excess,m_ipDemand,m_iVertNum*sizeof(Tflow));
	for (int v=0; v < m_iEdgeNum; ++v) {
		m_ipFlow[v]=delta=getLoCap(v);
		Excess[m_ipTail[v]]+=delta;
		Excess[m_ipHead[v]]-=delta;
	}
	m_ipPrec[m_iVertNum]=NIL;
	m_ipDepth[m_iVertNum]=0;
	m_ipThread[m_iVertNum]=0;
	m_ipPrice[m_iVertNum]=0;

	int edge=m_iEdgeNum;
	TRAN_setBigM();
	for (int v=0; v < m_iVertNum; ++v) {
		m_ipPrec[v]=edge;
		m_ipDepth[v]=1;
		m_ipThread[v]=v+1;
		if (isNonNegative(Excess[v])) {
			m_ipHead[edge]=v;
			m_ipTail[edge]=m_iVertNum;
			if (m_iProblemMsk & m_iLoCapMsk)
				m_ipLCap[edge]=0;
			m_ipUCap[edge]=m_ipFlow[edge]=Excess[v];
			m_ipCost[edge++]=m_iBigM;
			m_ipPrice[v]=m_iBigM;
		}
		else { // (Excess[v] < 0)
			m_ipHead[edge]=m_iVertNum;
			m_ipTail[edge]=v;
			if (m_iProblemMsk & m_iLoCapMsk)
				m_ipLCap[edge]=0;
			m_ipUCap[edge]=m_ipFlow[edge]=-Excess[v];
			m_ipCost[edge++]=m_iBigM;			
			m_ipPrice[v]=-m_iBigM;
		}		
	}

	m_iEdgeNum+=m_iVertNum;
	++m_iVertNum;
}	// end of CFlowNet<Tflow>::TRAN_initSolution

template <class Tflow>
void CFlowNet<Tflow>::TRAN_simplex(bool initSol)
{
	if (!initSol)
		TRAN_initSolution();

	int direct, side;
	int edge, DelEdge;

	for (; (edge=TRAN_checkFlowOpt(direct)) != NIL; ) {	// main loop
		if ((DelEdge= TRAN_changeFlow(edge,side,direct)) != edge) {
			TRAN_delEdge(DelEdge);
			TRAN_addEdge(edge,direct*side);
			TRAN_computePrices();
		}
	}

	m_iTRAN_Feasible=true;
	if (!initSol) {
		m_iVertNum--;
		m_iEdgeNum-=m_iVertNum;
		for (int e=0; e < m_iVertNum; ++e)
			if (isPositive(m_ipFlow[m_iEdgeNum+e])) {
				m_iTRAN_Feasible=false;
				break;
			}
	}

}	// end of CFlowNet<Tflow>::TRAN_simplex

template <class Tflow>
void CFlowNet<Tflow>::TRAN_printSolution()
{	
	printFlow();
	TRAN_tree();
}  // end of CFlowNet<Tflow>::PR_printSolution


template <class Tflow>
void CFlowNet<Tflow>::TRAN_tree()
{
	cout << "+-----------------------------------------+\n";
	cout << "|  Node  | Depth  |  Prec  | Thread |  Price |\n";
	for (int i=0; i < m_iVertNum; ++i) {
		cout << "| " << setw(6) << i
			 << " | " << setw(6) << m_ipDepth[i]
		     << " | " << setw(6) << m_ipPrec[i]
			<< " | " << setw(6) << m_ipThread[i]
			<< " | " << setw(6) << m_ipPrice[i] << "|\n";
	}
	cout << "+-----------------------------------------+\n";
}  // end of CFlowNet<Tflow>::TRAN_tree

template <class Tflow>
void CFlowNet<Tflow>::print()
{
	cout << "Flow:\n";
	for (int i=0; i < m_iEdgeNum-1; ++i)
		cout << m_ipFlow[i] << ", ";
	cout << m_ipFlow[m_iEdgeNum-1] << endl;

	cout << "Label:\n";
	for (int i=0; i < m_iVertNum-1; ++i)
		cout << m_ipLabel[i] << ", ";
	cout << m_ipLabel[m_iVertNum-1] << endl;

	cout << "Excess:\n";
	for (int i=0; i < m_iVertNum-1; ++i)
		cout << m_ipExcess[i] << ", ";
	cout << m_ipExcess[m_iVertNum-1] << endl;

	cout << "CurEdge:\n";
	for (int i=0; i < m_iVertNum-1; ++i)
		cout << m_ipCurEdge[i] << ", ";
	cout << m_ipCurEdge[m_iVertNum-1] << endl;

	cout << "Queue: First=" << m_iFirst << ", Last=" << m_iLast << endl;
	for (int i=0; i < m_iVertNum-1; ++i)
		cout << m_ipQueue[i] << ", ";
	cout << m_ipQueue[m_iVertNum-1] << endl;
}

template <class Tflow>
void CFlowNet<Tflow>::printNet()
{
	std::ofstream fout("Net.txt");
	fout << "Node: adjacency list (edge,other end\n";
	int e,v;
	for (int i=0; i < m_iVertNum; i++) {
		fout << i << ": ";
		for (int j=m_ipEdgeSep[i];j < m_ipEdgeSep[i+1]; j++) {
			e=m_ipEdge[j];
			v=(m_ipTail[e] == i)? m_ipHead[e]: m_ipTail[e];
			fout << "(" << e << "," << v << ")";
		}
		fout << "\n";
	}
	fout.close();
}	// end of CFlowNet<Tflow>::PrintNet

//////////////////////////////////////////////

template class CFlowNet<double>;
