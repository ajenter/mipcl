#include <cstring>
#include <cmath>
#include <iomanip>
#include <except.h>
#include "optportfolio.h"

Coptportfolio::Coptportfolio(const char* name): CMIP(name)
{
	char fileName[128];
	strcpy(fileName,name);
	strcat(fileName,".txt");
	readData(fileName);
}

#ifndef __ONE_THREAD_
Coptportfolio::Coptportfolio(const Coptportfolio &other, int thread): CMIP(other,thread)
{
} // end of Coptportfolio::Coptportfolio(Coptportfolio &other, int thread)

CMIP* Coptportfolio::clone(const CMIP *pMip, int thread)
{
	return static_cast<CMIP*>(new Coptportfolio(*static_cast<Coptportfolio*>(const_cast<CMIP*>(pMip)),thread));
}
#endif

Coptportfolio::~Coptportfolio()
{
#ifndef __ONE_THREAD_
	if (!m_iThread) {
#endif
		if (m_dpL)
			delete[] m_dpL;
#ifndef __ONE_THREAD_
	}
#endif
} // end of Coptportfolio::~Coptportfolio

//////////////////// implementation
void Coptportfolio::readData(const char* fileName)
{
	std::ifstream fin(fileName);
	if (!fin.is_open()) {
		throw new CFileException("Coptportfolio::readData",fileName);
	}

	int n, T, q;

	fin >> n >> m_iQ1 >> m_iQ2 >> T;
	m_iN=n; m_iT=T;
	if (!(m_dpL = new(std::nothrow) double[n*(T+2)+T])) {
		throw CMemoryException("Coptportfolio::readData");
	}
	m_dpP=(m_dpMu=m_dpL+n)+n;
	m_dpRet=m_dpP+T;
	
	for (int j=0; j < n; ++j)
		fin >> m_dpL[j];
	
	for (int t=0; t < T; ++t) {
		fin >> q;
		for (int j=0; j < n; ++j)
			fin >> m_dpRet[j*T+t];
	}
	fin.close();
} // end of Coptportfolio::readData()

void Coptportfolio::model(double p, double V)
{
	int r, n{m_iN}, T{m_iT}, q1{m_iQ1}, q2{m_iQ2};

	auto mu = [this](int j) {return this->m_dpMu[j];};
	auto l = [this](int j) {return this->m_dpL[j];};
	auto a = [this,T](int j, int t) {return sqrt(this->m_dpP[t])*this->m_dpRet[j*T+t];};

	auto x = [](int j) {return j;};
	auto v = [n](int j) {return j+n;};
	auto y = [n](int t) {return t+2*n;};

	computeParameters(p);
	m_dV = V;

	openMatrix(2*n+T+3,2*n+T+1,n*(7+T)+T,true,false);
	setObjSense(false);

	for (int j=0; j < n; ++j) { // add x-variables
		addVar(x(j),0,0.0,0.0,1.0);
	}
	for (int j=0; j < n; ++j) { // add v-variables
		addVar(v(j),VAR_BIN,0.0,0.0,1.0);
	}
	for (int t=0; t < T; ++t) { // add y-variables
		addVar(y(t),0,0.0,-VAR_INF,VAR_INF);
	}
	addVar(r=2*n+T,0,1.0,0.0,VAR_INF);

	addCtr(0,0,1.0,1.0); // sum{j=1}^n x_j = 1
	for (int j=0; j < n; ++j) {
		addEntry(1.0,0,x(j));
	}
	addCtr(1,0,V,INF); // sum{j=1}^n \mu_j x_j \ge ret
	for (int j=0; j < n; ++j) {
		addEntry(mu(j),1,x(j));
	}

	addCtr(2,0,q1,q2); // \sum_{j=1}^n v_j \qe q
	for (int j=0; j < n; ++j) {
		addEntry(1.0,2,v(j));
	}

	int row=2;
	for (int t=0; t < T; ++t) {
		m_ipArray[t]=y(t);
		addCtr(++row,0,0.0,0.0); // \sum_{j=1}^n a_{jt} x_j - y_t = 0
		for (int j=0; j < n; ++j) {
			addEntry(a(j,t),row,x(j));
		}
		addEntry(-1.0,row,y(t));
	}

	for (int j=0; j < n; ++j) {
		addCtr(++row,0,0.0,INF); // x_j - l_j v_j \ge 0
		addEntry(1.0,row,x(j));
		addEntry(-l(j),row,v(j));
		addCtr(++row,0,-INF,0.0); // x_j - v_j \le 0
		addEntry(1.0,row,x(j));
		addEntry(-1.0,row,v(j));
	}

	allowNormCtrs(1,T);

	addNormCtr(r,T,m_ipArray,0.0001);

	preprocOff();

	closeMatrix();

} // end of Coptportfolio::model()

void Coptportfolio::computeParameters(double p)
{
	double nom, den, f;
	double *dpR{m_dpRet}, *dpP{m_dpP};
	int n{m_iN}, T{m_iT};

	for (int j=0; j < n; ++j) {
		nom=den=0.0; f=1.0;
		for (int t=T; --t >= 0;) { // compute average rewards mu[j]
			nom+=f*log(dpR[t]);
			den+=f;
			f*=p;
		}
		m_dpMu[j]=exp(nom/den);
		dpR+=T;
	}
	den=0.0; f=1.0;
	for (int t=T; --t >= 0;) {
		dpP[t]=f;
		den+=f;
		f*=p;
	}
	for (int t=0; t < T; ++t) {
		dpP[t]/=den;
	}

// compute discounted revards ret(j,t)
	dpR=m_dpRet;
	for (int j=0; j < n; ++j) {
		den=0.0;
		for(int tau=0; tau < T; ++tau) {
			den+=dpP[tau]*dpR[tau];
		}
		for(int t=0; t < T; ++t) {
			dpR[t]-=den;
		}
		dpR+=T;
	}
} // end of Coptportfolio::computeParameters()

void Coptportfolio::printSolution(const char* fileName)
{
	if (getSolNum()) {
		std::ofstream fout(fileName);
		if (!fout.is_open()) {
			throw new CFileException("Coptportfolio::printSolution()",fileName);
		}

		fout.precision(4);
		fout << "Revenue: " << m_dV << ", Risk: " << getObjVal() << std::endl;
		for (int j=0; j <= m_iN; ++j) {
			fout << "=========";
		}
		fout << std::endl;
		fout.setf(std::ios_base::right,std::ios_base::adjustfield);
		for (int j=1; j <= m_iN; ++j) {
			fout << std::setw(10) << j;
		}
		fout << std::endl;
		for (int j=0; j <= m_iN; ++j) {
			fout << "=========";
		}
		fout << std::endl;
		fout.setf(std::ios_base::fixed,std::ios_base::floatfield);
		double *dpX;
		int *ipHd;
		int n=getSolution(dpX,ipHd);

		for (int j=0; j < m_iN; ++j) {
			fout << std::setw(10) << dpX[j];
		}
		fout << std::endl;
		fout.close();
	}
} // end of Coptportfolio::printSolution()

