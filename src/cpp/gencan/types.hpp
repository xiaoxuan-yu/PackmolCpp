#pragma once

enum class GencanImplMode {
    kFortran = 0,
    kCpp = 1,
};

struct GencanProblemView {
    const int* n;
    double* x;
    const double* l;
    const double* u;
    const int* m;
    const double* lambda;
    const double* rho;
    const int* gtype;
    const int* htvtype;
    const int* trtype;
    const int* iprint;
    const int* ncomp;
};

struct GencanOptions {
    const double* epsgpen;
    const double* epsgpsn;
    const int* maxitnfp;
    const double* epsnfp;
    const int* maxitngp;
    const double* fmin;
    const int* maxit;
    const int* maxfc;
    const double* udelta0;
    const int* ucgmaxit;
    const int* cgscre;
    const double* cggpnf;
    const double* cgepsi;
    const double* cgepsf;
    const double* epsnqmp;
    const int* maxitnqmp;
    const bool* nearlyq;
    const double* nint;
    const double* next;
    const int* mininterp;
    const int* maxextrap;
    const double* eta;
    const double* delmin;
    const double* lspgma;
    const double* lspgmi;
    const double* theta;
    const double* gamma;
    const double* beta;
    const double* sigma1;
    const double* sigma2;
    const double* sterel;
    const double* steabs;
    const double* epsrel;
    const double* epsabs;
    const double* infrel;
    const double* infabs;
};

struct GencanCounters {
    int* iter;
    int* fcnt;
    int* gcnt;
    int* cgcnt;
    int* spgiter;
    int* spgfcnt;
    int* tniter;
    int* tnfcnt;
    int* tnstpcnt;
    int* tnintcnt;
    int* tnexgcnt;
    int* tnexbcnt;
    int* tnintfe;
    int* tnexgfe;
    int* tnexbfe;
};

struct GencanMutableState {
    double* f;
    double* g;
    double* gpeucn2;
    double* gpsupn;
    int* inform;
    int* ind;
    double* lastgpns;
};

struct GencanWorkspace {
    double* s;
    double* y;
    double* d;
    double* w;
};

struct GencanSpglsContext {
    const int* n;
    double* x;
    const int* m;
    const double* lambda;
    const double* rho;
    double* f;
    const double* g;
    const double* l;
    const double* u;
    const double* lamspg;
    const double* nint;
    const int* mininterp;
    const double* fmin;
    const int* maxfc;
    const int* iprint;
    int* fcnt;
    int* inform;
    double* xtrial;
    double* d;
    const double* gamma;
    const double* sigma1;
    const double* sigma2;
    const double* sterel;
    const double* steabs;
    const double* epsrel;
    const double* epsabs;
    const double* infrel;
    const double* infabs;
};

struct GencanTnlsContext {
    const int* nind;
    const int* ind;
    const int* n;
    double* x;
    const int* m;
    const double* lambda;
    const double* rho;
    const double* l;
    const double* u;
    double* f;
    double* g;
    const double* d;
    const double* amax;
    const int* rbdtype;
    const int* rbdind;
    const double* nint;
    const double* next;
    const int* mininterp;
    const int* maxextrap;
    const double* fmin;
    const int* maxfc;
    const int* gtype;
    const int* iprint;
    int* fcnt;
    int* gcnt;
    int* intcnt;
    int* exgcnt;
    int* exbcnt;
    int* inform;
    double* xplus;
    double* xtmp;
    double* xbext;
    const double* gamma;
    const double* beta;
    const double* sigma1;
    const double* sigma2;
    const double* sterel;
    const double* steabs;
    const double* epsrel;
    const double* epsabs;
    const double* infrel;
    const double* infabs;
};

struct GencanCgContext {
    const int* nind;
    const int* ind;
    const int* n;
    double* x;
    const int* m;
    const double* lambda;
    const double* rho;
    double* g;
    const double* delta;
    const double* l;
    const double* u;
    const double* eps;
    const double* epsnqmp;
    const int* maxitnqmp;
    const int* maxit;
    const bool* nearlyq;
    const int* gtype;
    const int* htvtype;
    const int* trtype;
    const int* iprint;
    const int* ncomp;
    double* s;
    int* iter;
    int* rbdtype;
    int* rbdind;
    int* inform;
    double* w;
    double* y;
    double* r;
    double* d;
    double* sprev;
    const double* theta;
    const double* sterel;
    const double* steabs;
    const double* epsrel;
    const double* epsabs;
    const double* infrel;
    const double* infabs;
};
