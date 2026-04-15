#pragma once

#include "gencan/types.hpp"

#include <vector>

extern thread_local int g_evalhd_fortran_call_count;

struct PackmolPrecisionState {
    double precision_value;
    double fdist_value;
    double frest_value;
};

struct PackmolEvalMetaState {
    int maxrest;
    int mrperatom;
    int ntype;
    int ntotmol;
    int ntotat;
    int nfixedat;
    int ntype_with_fixed;
    int lcellfirst;
    int ncells[3];
    double scale;
    double scale2;
    double pbc_min[3];
    double pbc_length[3];
    double cell_length[3];
    bool init1;
    bool move;
    bool fix;
    bool using_pbc;
};

struct PackmolEvalRawPointers {
    const int* nmols;
    const int* natoms;
    const int* idfirst;
    const int* nratom;
    const int* iratom;
    const int* ityperest;
    int* ibmol;
    int* ibtype;
    int* latomnext;
    int* latomfirst;
    int* lcellnext;
    double* xcart;
    const double* coor;
    const double* restpars;
    const double* radius;
    const double* radius_ini;
    const double* fscale;
    const double* short_radius;
    const double* short_radius_scale;
    double* gxcar;
    double* fdist_atom;
    double* frest_atom;
};

struct PackmolEvalStateSnapshot {
    PackmolEvalMetaState meta;
    PackmolEvalRawPointers raw;
    std::vector<int> comptype_flags;
    std::vector<int> fixedatom_flags;
    std::vector<int> use_short_radius_flags;
};

GencanImplMode active_impl_mode();
bool gencan_debug_enabled();
bool cg_dtw_relax_enabled();
bool use_cpp_numeric_kernel();
bool use_cpp_eval_kernel();
void run_spgls_context_cpp(const GencanSpglsContext& context);
void run_tnls_context_cpp(const GencanTnlsContext& context);
void run_cg_context_cpp(const GencanCgContext& context);
double norm2_kernel(int n, const double* x);

int packmol_gencan_cpp_probe_impl(int n, const double* x, double* fx_out);
int packmol_gencan_impl_mode_cpp();
void packmol_gencan_reset_evalhd_call_count_cpp();
int packmol_gencan_evalhd_call_count_cpp();

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

void eval_objective_full_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
);

bool try_eval_objective_full_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
);

bool eval_objective_full_cpp_direct(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
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

PackmolPrecisionState read_packmol_precision_state_cpp();
bool load_packmol_eval_state_cpp(PackmolEvalStateSnapshot* state);
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
);

void packmol_gencan_gencan_solver_cpp(
    const int* n,
    double* x,
    const double* l,
    const double* u,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* epsgpen,
    const double* epsgpsn,
    const int* maxitnfp,
    const double* epsnfp,
    const int* maxitngp,
    const double* fmin,
    const int* maxit,
    const int* maxfc,
    const double* udelta0,
    const int* ucgmaxit,
    const int* cgscre,
    const double* cggpnf,
    const double* cgepsi,
    const double* cgepsf,
    const double* epsnqmp,
    const int* maxitnqmp,
    const bool* nearlyq,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* f,
    double* g,
    double* gpeucn2,
    double* gpsupn,
    int* iter,
    int* fcnt,
    int* gcnt,
    int* cgcnt,
    int* spgiter,
    int* spgfcnt,
    int* tniter,
    int* tnfcnt,
    int* tnstpcnt,
    int* tnintcnt,
    int* tnexgcnt,
    int* tnexbcnt,
    int* tnintfe,
    int* tnexgfe,
    int* tnexbfe,
    int* inform,
    double* s,
    double* y,
    double* d,
    int* ind,
    double* lastgpns,
    double* w,
    const double* eta,
    const double* delmin,
    const double* lspgma,
    const double* lspgmi,
    const double* theta,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);
