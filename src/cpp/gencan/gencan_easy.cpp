#include "gencan/api.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

extern "C" void packmol_easyg_fortran_c(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpsn,
    const int* maxit,
    const int* maxfc,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* inform,
    int* wi,
    double* wd,
    double* delmin
);

namespace {

void easygencan_cpp(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpsn,
    const int* maxit,
    const int* maxfc,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* inform,
    int* wi,
    double* wd,
    double* delmin
) {
    const int n_val = *n;
    const int tmax = 1000;

    const double steabs = 1.0e-10;
    const double sterel = 1.0e-7;
    const double epsabs = 1.0e-20;
    const double epsrel = 1.0e-10;
    const double infabs = 1.0e99;
    const double infrel = 1.0e20;

    const double beta = 0.5;
    const double gamma = 1.0e-4;
    const double theta = 1.0e-6;
    const double sigma1 = 0.1;
    const double sigma2 = 0.9;

    const int maxextrap = 100;
    const int mininterp = 4;
    const double nint = 2.0;
    const double next = 2.0;
    double delmin_local = 1.0e-2;
    const double eta = 0.9;
    const double lspgma = 1.0e10;
    const double lspgmi = 1.0e-10;

    const int gtype = 0;
    const int htvtype = 1;
    const double epsgpen = 0.0;

    const int maxitngp = tmax;
    const int maxitnfp = *maxit;
    const double epsnfp = 0.0;
    const double fmin = 1.0e-5;

    const double delta0 = -1.0;
    const int cgmaxit = -1;
    const int trtype = 1;
    const int cgscre = 2;
    const double cggpnf = std::max(1.0e-4, std::max(epsgpen, *epsgpsn));
    const double cgepsi = 1.0e-1;
    const double cgepsf = 1.0e-5;
    const double epsnqmp = 1.0e-4;
    const int maxitnqmp = 5;
    const bool nearlyq = false;

    int spgiter = 0;
    int spgfcnt = 0;
    int tniter = 0;
    int tnfcnt = 0;
    int tnstpcnt = 0;
    int tnintcnt = 0;
    int tnexgcnt = 0;
    int tnexbcnt = 0;
    int tnintfe = 0;
    int tnexgfe = 0;
    int tnexbfe = 0;
    double gpeucn2 = 0.0;

    std::vector<double> lastgpns(tmax, 0.0);

    double* s = wd;
    double* y = wd + n_val;
    double* d = wd + 2 * n_val;
    double* w = wd + 3 * n_val;

    packmol_gencan_gencan_solver_cpp(
        n, x, l, u, m, lambda, rho, &epsgpen, epsgpsn, &maxitnfp, &epsnfp,
        &maxitngp, &fmin, maxit, maxfc, &delta0, &cgmaxit, &cgscre, &cggpnf,
        &cgepsi, &cgepsf, &epsnqmp, &maxitnqmp, &nearlyq, &nint, &next,
        &mininterp, &maxextrap, &gtype, &htvtype, &trtype, iprint, ncomp, f, g,
        &gpeucn2, gpsupn, iter, fcnt, gcnt, cgcnt, &spgiter, &spgfcnt, &tniter,
        &tnfcnt, &tnstpcnt, &tnintcnt, &tnexgcnt, &tnexbcnt, &tnintfe, &tnexgfe,
        &tnexbfe, inform, s, y, d, wi, lastgpns.data(), w, &eta, &delmin_local, &lspgma,
        &lspgmi, &theta, &gamma, &beta, &sigma1, &sigma2, &sterel, &steabs,
        &epsrel, &epsabs, &infrel, &infabs
    );

    if (delmin != nullptr) {
        *delmin = delmin_local;
    }
}

} // namespace

void packmol_gencan_easy_bridge_impl(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpsn,
    const int* maxit,
    const int* maxfc,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* inform,
    int* wi,
    double* wd,
    double* delmin
) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_easyg_fortran_c(
                n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, trtype, iprint,
                ncomp, f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin
            );
            return;
        case GencanImplMode::kCpp:
            easygencan_cpp(
                n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, iprint, ncomp,
                f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin
            );
            return;
    }
}
