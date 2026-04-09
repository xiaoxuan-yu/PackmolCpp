#pragma once

extern thread_local int g_evalhd_fortran_call_count;

void shrink_inplace(int nind, const int* ind, double* v);
void expand_inplace(int nind, const int* ind, double* v);

void gp_ieee_signal1_cpp(
    double gpsupn,
    double* acgeps,
    double* bcgeps,
    double cgepsf,
    double cgepsi,
    double cggpnf
);

void gp_ieee_signal2_cpp(
    int* cgmaxit,
    int nind,
    bool nearlyq,
    int ucgmaxit,
    int cgscre,
    double* kappa,
    double gpeucn2,
    double gpeucn20,
    double epsgpen2,
    double epsgpsn,
    double* cgeps,
    double acgeps,
    double bcgeps,
    double cgepsf,
    double cgepsi,
    double gpsupn,
    double gpsupn0
);

int evaluate_post_step_inform_cpp(
    bool precision_after,
    int line_inform,
    double gpeucn2_after,
    double epsgpen,
    double gpsupn_after,
    double epsgpsn,
    double f_before,
    double f_after,
    double epsnfp,
    int maxitnfp,
    double infabs,
    const double* lastgpns,
    int maxitngp,
    double fmin,
    int iter_value,
    int maxit,
    int fcnt_value,
    int maxfc,
    double* progress_fprev,
    double* progress_bestprog,
    int* progress_itnfp
);

int finalize_post_step_inform_cpp(
    int post_inform,
    bool precision_after,
    double gpeucn2_after,
    double epsgpen,
    double gpsupn_after,
    double epsgpsn,
    int iter_value,
    int maxit,
    int fcnt_value,
    int maxfc
);

void calcf_cpp_reduced(
    const int* nind,
    const int* ind,
    double* x,
    const int* n,
    const double* xc,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
);

void calcg_cpp_reduced(
    const int* nind,
    const int* ind,
    double* x,
    const int* n,
    const double* xc,
    const int* m,
    const double* lambda,
    const double* rho,
    const int* gtype,
    double* g,
    const double* sterel,
    const double* steabs,
    int* inform
);

void eval_gradient_full_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const int* gtype,
    double* g,
    const double* sterel,
    const double* steabs,
    int* inform
);

bool packmolprecision_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    int* eval_flag
);

void calchddiff_cpp_reduced(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    double* d,
    double* g,
    const int* m,
    const double* lambda,
    const double* rho,
    const int* gtype,
    double* hd,
    double* xtmp,
    const double* sterel,
    const double* steabs,
    int* inform
);

void calchd_cpp_reduced(
    const int* nind,
    const int* ind,
    double* x,
    double* d,
    double* g,
    const int* n,
    const double* xc,
    double* hd
);
