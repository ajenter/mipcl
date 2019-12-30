#include <fstream>
#include <iostream>
#include <cstring>
#include <except.h>
#include "oneMachine.h"

int read(const char* fileName,
         int &n, int* &r, int* &d, int* &p, double* &c)
{
    std::ifstream fin(fileName);
    if (!fin.is_open()) {
        throw CFileException("read",fileName);
    }
    fin.ignore(256,'\n');
    fin >> n;
    r = new int[n];
    d = new int[n];
    p = new int[n];
    c = new double[n]; 
    for (int i=0; i < n; i++) {
        fin >> r[i] >> d[i] >> p[i] >> c[i];
    }
    return 0;
}

int main(int argc, char* argv[])
{
    int n, *r,*d,*p;
    double *c;
    if (argc < 2) {
        std::cerr << "Enter file name as the argument!\n"; 
        return 1;
    }
    try {
        read(argv[1],n,r,d,p,c);
        std::cerr << "n=" << n << std::endl;
        COneMachine prob("OneMachine",n,c,p,r,d);
        prob.optimize();
        char str[256];
        strcpy(str,argv[1]);
        strcat(str,".sol");
        prob.printSchedule(str);
    }
    catch(CException* pe) {
        std::cerr << pe->what() << std::endl;
        delete pe;
        return 2;
    }
    if (c) delete[] c;
    if (p) delete[] p;
    if (d) delete[] d;
    if (r) delete[] r;
    return 0;
}
