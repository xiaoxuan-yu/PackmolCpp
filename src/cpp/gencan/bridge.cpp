#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

thread_local int g_cg_active_call_id = 0;
thread_local int g_spg_post_cpp_replay_depth = 0;
thread_local int g_tn_post_cpp_replay_depth = 0;
thread_local int g_init1_phase_active = 0;

enum class GencanImplMode {
    kFortran = 0,
    kCpp = 1,
    kAb = 2
};

GencanImplMode parse_impl_mode() {
    const char* env = std::getenv("PACKMOL_GENCAN_IMPL");
    if (env == nullptr) {
        return GencanImplMode::kCpp;
    }

    std::string value(env);
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }

    if (value == "fortran") {
        return GencanImplMode::kFortran;
    }
    if (value == "ab") {
        return GencanImplMode::kAb;
    }
    return GencanImplMode::kCpp;
}

GencanImplMode active_impl_mode() {
    static const GencanImplMode mode = parse_impl_mode();
    return mode;
}

bool gencan_debug_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_DEBUG");
    if (env == nullptr) {
        return false;
    }
    std::string value(env);
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value == "1" || value == "true" || value == "on" || value == "yes";
}

bool cg_shadow_compare_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_CG_SHADOW");
    if (env == nullptr) {
        return false;
    }
    std::string value(env);
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value == "1" || value == "true" || value == "on" || value == "yes";
}

bool cg_dtw_relax_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_CG_DTW_RELAX");
    if (env == nullptr) {
        return true;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool use_cpp_numeric_kernel() {
    const char* env = std::getenv("PACKMOL_GENCAN_NUMERIC_CPP");
    if (env == nullptr) {
        return true;
    }
    return !(env[0] == '0' || env[0] == 'f' || env[0] == 'F' || env[0] == 'n' || env[0] == 'N');
}

bool tn_post_shadow_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_SHADOW");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool fallback_seed_state_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_FALLBACK_SEED_STATE");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool spg_post_handoff_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_SPG_POST_CPP_HANDOFF");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool tn_post_handoff_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_CPP_HANDOFF");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool tn_post_handoff_safe_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_SAFE");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool tn_post_handoff_unsafe_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool tn_post_handoff_canonicalize_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CANONICALIZE");
    if (env == nullptr) {
        return true;
    }
    return !(env[0] == '0' || env[0] == 'f' || env[0] == 'F' || env[0] == 'n' || env[0] == 'N');
}

bool tn_post_handoff_cpp_replay_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY");
    if (env == nullptr) {
        return true;
    }
    return !(env[0] == '0' || env[0] == 'f' || env[0] == 'F' || env[0] == 'n' || env[0] == 'N');
}

bool block_cpp_tail_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_BLOCK_CPP_TAIL");
    if (env == nullptr) {
        // Legacy alias kept only for transition; blocked-tail is fully C++ now.
        env = std::getenv("PACKMOL_GENCAN_BLOCK_FORTRAN_TAIL");
    }
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool cpp_tail_reduction_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_CPP_TAIL_REDUCTION");
    if (env == nullptr) {
        return true;
    }
    return !(env[0] == '0' || env[0] == 'f' || env[0] == 'F' || env[0] == 'n' || env[0] == 'N');
}

bool ab_compare_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_AB_COMPARE");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool block_cpp_tail_enabled();
bool init1_phase_active_cpp();

bool handoff_cpp_tail_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_HANDOFF_CPP_TAIL");
    if (env == nullptr) {
        return (active_impl_mode() == GencanImplMode::kCpp && !ab_compare_enabled()) ||
               block_cpp_tail_enabled();
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool handoff_cpp_tail_allow_init1_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_HANDOFF_CPP_TAIL_ALLOW_INIT1");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool handoff_cpp_tail_allow_spg_post_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_HANDOFF_CPP_TAIL_ALLOW_SPG_POST");
    if (env == nullptr) {
        return active_impl_mode() == GencanImplMode::kCpp && !ab_compare_enabled();
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool handoff_cpp_tail_allow_tn_post_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_HANDOFF_CPP_TAIL_ALLOW_TN_POST");
    if (env == nullptr) {
        return active_impl_mode() == GencanImplMode::kCpp && !ab_compare_enabled();
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool handoff_cpp_tail_allow_tn_post_explicitly_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_HANDOFF_CPP_TAIL_ALLOW_TN_POST");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool is_tn_tail_reason_cpp(const std::string& reason) {
    return reason == "tn_post_nonterminal" || reason == "tn_no_free_variables";
}

int spg_post_cpp_replay_max_depth() {
    const char* env = std::getenv("PACKMOL_GENCAN_SPG_POST_CPP_REPLAY_MAX_DEPTH");
    if (env == nullptr || env[0] == '\0') {
        return 1;
    }
    char* end = nullptr;
    long value = std::strtol(env, &end, 10);
    if (end == env || value <= 0) {
        return 1;
    }
    if (value > 4) {
        return 4;
    }
    return static_cast<int>(value);
}

int tn_post_cpp_replay_max_depth() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_CPP_REPLAY_MAX_DEPTH");
    if (env == nullptr || env[0] == '\0') {
        return 1;
    }
    char* end = nullptr;
    long value = std::strtol(env, &end, 10);
    if (end == env || value <= 0) {
        return 1;
    }
    if (value > 256) {
        return 256;
    }
    return static_cast<int>(value);
}

bool tn_post_cpp_replay_tail_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_CPP_REPLAY_TAIL");
    if (env == nullptr) {
        return active_impl_mode() == GencanImplMode::kCpp && !ab_compare_enabled();
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

void apply_spg_post_cpp_replay_bookkeeping(
    int* iter,
    int* fcnt,
    int* gcnt
) {
    // Replay starts from the post-line-search state and skips one outer
    // accounting window present in the legacy continuation path.
    *iter += 1;
    *fcnt += 2;
    *gcnt += 1;
}

int block_tail_continuation_steps() {
    const char* env = std::getenv("PACKMOL_GENCAN_BLOCK_TAIL_STEPS");
    if (env == nullptr || env[0] == '\0') {
        // Default: keep continuing until a regular GENCAN stop condition
        // is reached (maxit/maxfc/post-inform), matching main-loop intent.
        return 0;
    }
    char* end = nullptr;
    long value = std::strtol(env, &end, 10);
    if (end == env || value <= 0) {
        return 0;
    }
    if (value > 1024) {
        return 1024;
    }
    return static_cast<int>(value);
}

int tn_post_retry_spg_steps() {
    const char* env = std::getenv("PACKMOL_GENCAN_TN_POST_RETRY_SPG");
    if (env == nullptr || env[0] == '\0') {
        return 0;
    }
    char* end = nullptr;
    long value = std::strtol(env, &end, 10);
    if (end == env || value <= 0) {
        return 0;
    }
    if (value > 16) {
        return 16;
    }
    return static_cast<int>(value);
}

int spg_post_retry_steps() {
    const char* env = std::getenv("PACKMOL_GENCAN_SPG_POST_RETRY");
    if (env == nullptr || env[0] == '\0') {
        return 0;
    }
    char* end = nullptr;
    long value = std::strtol(env, &end, 10);
    if (end == env || value <= 0) {
        return 0;
    }
    if (value > 16) {
        return 16;
    }
    return static_cast<int>(value);
}

void shrink_inplace(const int nind, const int* ind, double* v) {
    for (int i = 1; i <= nind; ++i) {
        const int indi = ind[i - 1];
        if (i != indi) {
            const double tmp = v[indi - 1];
            v[indi - 1] = v[i - 1];
            v[i - 1] = tmp;
        }
    }
}

void expand_inplace(const int nind, const int* ind, double* v) {
    for (int i = nind; i >= 1; --i) {
        const int indi = ind[i - 1];
        if (i != indi) {
            const double tmp = v[indi - 1];
            v[indi - 1] = v[i - 1];
            v[i - 1] = tmp;
        }
    }
}

void gp_ieee_signal1_cpp(
    const double gpsupn,
    double* acgeps,
    double* bcgeps,
    const double cgepsf,
    const double cgepsi,
    const double cggpnf
) {
    if (gpsupn > 0.0) {
        *acgeps = std::log10(cgepsf / cgepsi) / std::log10(cggpnf / gpsupn);
        *bcgeps = std::log10(cgepsi) - (*acgeps) * std::log10(gpsupn);
    } else {
        *acgeps = 0.0;
        *bcgeps = cgepsf;
    }
}

void gp_ieee_signal2_cpp(
    int* cgmaxit,
    const int nind,
    const bool nearlyq,
    const int ucgmaxit,
    const int cgscre,
    double* kappa,
    const double gpeucn2,
    const double gpeucn20,
    const double epsgpen2,
    const double epsgpsn,
    double* cgeps,
    const double acgeps,
    const double bcgeps,
    const double cgepsf,
    const double cgepsi,
    const double gpsupn,
    const double gpsupn0
) {
    if (ucgmaxit <= 0) {
        if (nearlyq) {
            *cgmaxit = nind;
        } else {
            if (cgscre == 1) {
                *kappa = std::log10(gpeucn2 / gpeucn20) / std::log10(epsgpen2 / gpeucn20);
            } else {
                *kappa = std::log10(gpsupn / gpsupn0) / std::log10(epsgpsn / gpsupn0);
            }
            *kappa = std::max(0.0, std::min(1.0, *kappa));
            *cgmaxit = static_cast<int>(
                (1.0 - *kappa) * std::max(1.0, 10.0 * std::log10(static_cast<double>(nind))) +
                (*kappa) * static_cast<double>(nind)
            );
            *cgmaxit = std::min(20, *cgmaxit);
        }
    } else {
        *cgmaxit = ucgmaxit;
    }

    if (cgscre == 1) {
        *cgeps = std::sqrt(std::pow(10.0, acgeps * std::log10(gpeucn2) + bcgeps));
    } else {
        *cgeps = std::pow(10.0, acgeps * std::log10(gpsupn) + bcgeps);
    }
    *cgeps = std::max(cgepsf, std::min(cgepsi, *cgeps));
}

int evaluate_post_step_inform_cpp(
    const bool precision_after,
    const int line_inform,
    const double gpeucn2_after,
    const double epsgpen,
    const double gpsupn_after,
    const double epsgpsn,
    const double f_before,
    const double f_after,
    const double epsnfp,
    const int maxitnfp,
    const double infabs,
    const double* lastgpns,
    const int maxitngp,
    const double fmin,
    const int iter_value,
    const int maxit,
    const int fcnt_value,
    const int maxfc,
    double* progress_fprev,
    double* progress_bestprog,
    int* progress_itnfp
) {
    int post_inform = -1;
    if (precision_after) {
        return line_inform;
    }

    if (gpeucn2_after <= epsgpen * epsgpen) {
        post_inform = 0;
    } else if (gpsupn_after <= epsgpsn) {
        post_inform = 1;
    } else {
        double currprog = f_before - f_after;
        if (progress_fprev != nullptr) {
            currprog = (*progress_fprev) - f_after;
        }
        double bestprog = std::max(currprog, 0.0);
        if (progress_bestprog != nullptr) {
            *progress_bestprog = std::max(currprog, *progress_bestprog);
            bestprog = *progress_bestprog;
        }
        int itnfp = 0;
        if (progress_itnfp != nullptr) {
            itnfp = *progress_itnfp;
        }
        if (currprog <= epsnfp * bestprog) {
            itnfp += 1;
            if (itnfp >= maxitnfp) {
                post_inform = 2;
            }
        } else {
            itnfp = 0;
        }
        if (progress_itnfp != nullptr) {
            *progress_itnfp = itnfp;
        }
        if (progress_fprev != nullptr) {
            *progress_fprev = f_after;
        }
    }

    if (post_inform < 0) {
        double gpnmax = infabs;
        for (int i = 0; i < maxitngp; ++i) {
            gpnmax = std::max(gpnmax, lastgpns[i]);
        }
        if (gpeucn2_after >= gpnmax) {
            post_inform = 3;
        }
    }

    if (post_inform < 0 && f_after <= fmin) {
        post_inform = 4;
    } else if (post_inform < 0 && iter_value >= maxit) {
        post_inform = 7;
    } else if (post_inform < 0 && fcnt_value >= maxfc) {
        post_inform = 8;
    }

    if (post_inform < 0 && line_inform == 6) {
        post_inform = 6;
    }

    return post_inform;
}

void project_state_and_metrics_cpp(
    const int n_val,
    double* x,
    const double* g,
    const double* l,
    const double* u,
    const double epsrel,
    const double epsabs,
    int* ind,
    int* nind_after,
    double* gpsupn_after,
    double* gpeucn2_after
) {
    *gpsupn_after = 0.0;
    *gpeucn2_after = 0.0;
    *nind_after = 0;
    for (int i = 0; i < n_val; ++i) {
        if (x[i] <= l[i] + std::max(epsrel * std::abs(l[i]), epsabs)) {
            x[i] = l[i];
        } else if (x[i] >= u[i] - std::max(epsrel * std::abs(u[i]), epsabs)) {
            x[i] = u[i];
        }
        const double xpg = x[i] - g[i];
        const double gpi = std::min(u[i], std::max(l[i], xpg)) - x[i];
        *gpsupn_after = std::max(*gpsupn_after, std::abs(gpi));
        *gpeucn2_after += gpi * gpi;
        if (x[i] > l[i] && x[i] < u[i]) {
            ind[*nind_after] = i + 1;
            *nind_after += 1;
        }
    }
}

int finalize_blocked_tail_inform_cpp(
    const int post_inform,
    const bool precision_after,
    const double gpeucn2_after,
    const double epsgpen,
    const double gpsupn_after,
    const double epsgpsn,
    const int iter_value,
    const int maxit,
    const int fcnt_value,
    const int maxfc
) {
    if (post_inform >= 0) {
        return post_inform;
    }
    if (precision_after || gpeucn2_after <= epsgpen * epsgpen) {
        return 0;
    }
    if (gpsupn_after <= epsgpsn) {
        return 1;
    }
    if (iter_value >= maxit) {
        return 7;
    }
    if (fcnt_value >= maxfc) {
        return 8;
    }
    return 6;
}

double norm2_kernel(const int n, const double* x);

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

void spgls_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    const double* g,
    const double* l,
    const double* u,
    const double* lamspg,
    const double* nint,
    const int* mininterp,
    const double* fmin,
    const int* maxfc,
    int* fcnt,
    int* inform,
    double* xtrial,
    double* d,
    const double* gamma,
    const double* sigma1,
    const double* sigma2,
    const double* epsrel,
    const double* epsabs
);

bool tnls_cpp_subset(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* l,
    const double* u,
    double* f,
    double* g,
    const double* d,
    const double* amax,
    const int* rbdtype,
    const int* rbdind,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    const int* maxfc,
    const int* gtype,
    int* fcnt,
    int* gcnt,
    int* intcnt,
    int* exgcnt,
    int* exbcnt,
    int* inform,
    double* xplus,
    double* xtmp,
    double* xbext,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs
);

void cg_cpp_full(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* eps,
    const double* epsnqmp,
    const int* maxitnqmp,
    const int* maxit,
    const bool* nearlyq,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const double* theta,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infabs,
    double* s,
    double* w,
    double* y,
    double* r,
    double* d,
    double* sprev,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform
);

extern "C" void packmol_gencan_cg_bridge(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* eps,
    const double* epsnqmp,
    const int* maxitnqmp,
    const int* maxit,
    const bool* nearlyq,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* s,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform,
    double* w,
    double* y,
    double* r,
    double* d,
    double* sprev,
    const double* theta,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);

struct TnPostContinuationInputs {
    const int* n;
    double* x;
    const int* m;
    const double* lambda;
    const double* rho;
    const int* gtype;
    double* g;
    const double* l;
    const double* u;
    const double* epsrel;
    const double* epsabs;
    int* ind;
    double* gpeucn2;
    double* gpsupn;
    const int* maxitngp;
    double* lastgpns;
    int* iter;
    const int* maxfc;
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
    double* s_vec;
    double* y_vec;
    const double* udelta0;
    const int* ucgmaxit;
    const int* cgscre;
    const double* cggpnf;
    const double* cgepsi;
    const double* cgepsf;
    const double* epsnqmp;
    const int* maxitnqmp;
    const bool* nearlyq;
    const int* htvtype;
    const int* trtype;
    const int* iprint;
    const int* ncomp;
    const double* eta;
    const double* delmin;
    const double* lspgma;
    const double* lspgmi;
    const double* nint;
    const double* next;
    const int* mininterp;
    const int* maxextrap;
    const double* fmin;
    double* f;
    const double* theta;
    const double* gamma;
    const double* beta;
    const double* sigma1;
    const double* sigma2;
    const double* sterel;
    const double* steabs;
    const double* epsgpen;
    const double* epsgpsn;
    const int* maxitnfp;
    const double* epsnfp;
    const double* infrel;
    const double* infabs;
    const int* maxit;
    const bool cg_schedule_seed_valid;
    const double cg_schedule_acgeps;
    const double cg_schedule_bcgeps;
    const double cg_schedule_gpeucn20;
    const double cg_schedule_gpsupn0;
    int* inform;
};

void continue_tn_post_blocked_tail_cpp(
    const int n_val,
    const TnPostContinuationInputs& inputs,
    bool* precision_after,
    double* gpsupn_after,
    double* gpeucn2_after,
    int* nind_after,
    int* post_inform,
    int* continuation_budget_used,
    int* continuation_steps_done,
    double* progress_fprev,
    double* progress_bestprog,
    int* progress_itnfp
) {
    const int continuation_budget = block_tail_continuation_steps();
    const bool unlimited_continuation = continuation_budget <= 0;
    *continuation_budget_used = continuation_budget;
    std::vector<double> xtrial_work(n_val, 0.0);
    std::vector<double> d_work(n_val, 0.0);
    std::vector<double> x_prev(n_val, 0.0);
    std::vector<double> g_prev(n_val, 0.0);
    double sts = 0.0;
    double sty = 0.0;
    if (inputs.s_vec != nullptr && inputs.y_vec != nullptr) {
        for (int i = 0; i < n_val; ++i) {
            sts += inputs.s_vec[i] * inputs.s_vec[i];
            sty += inputs.s_vec[i] * inputs.y_vec[i];
        }
    }
    while (unlimited_continuation || *continuation_steps_done < continuation_budget) {
        *continuation_steps_done += 1;
        *inputs.iter += 1;
        const double f_before_retry = *inputs.f;
        for (int i = 0; i < n_val; ++i) {
            x_prev[i] = inputs.x[i];
            g_prev[i] = inputs.g[i];
        }
        int line_inform = 0;
        double gieucn2_after = 0.0;
        for (int i = 0; i < n_val; ++i) {
            const double xpg = inputs.x[i] - inputs.g[i];
            const double gpi = std::min(inputs.u[i], std::max(inputs.l[i], xpg)) - inputs.x[i];
            if (inputs.x[i] > inputs.l[i] && inputs.x[i] < inputs.u[i]) {
                gieucn2_after += gpi * gpi;
            }
        }
        const double ometa2 = (1.0 - *inputs.eta) * (1.0 - *inputs.eta);
        const bool use_spg_step =
            (*gpeucn2_after > 0.0 && gieucn2_after <= ometa2 * (*gpeucn2_after)) ||
            *nind_after <= 0;

        if (!use_spg_step) {
            std::vector<double> x_work(inputs.x, inputs.x + n_val);
            std::vector<double> g_work(inputs.g, inputs.g + n_val);
            std::vector<double> l_work(inputs.l, inputs.l + n_val);
            std::vector<double> u_work(inputs.u, inputs.u + n_val);
            const int nind_work = *nind_after;

            shrink_inplace(nind_work, inputs.ind, x_work.data());
            shrink_inplace(nind_work, inputs.ind, g_work.data());
            shrink_inplace(nind_work, inputs.ind, l_work.data());
            shrink_inplace(nind_work, inputs.ind, u_work.data());

            const double xnorm_full = std::sqrt(norm2_kernel(n_val, inputs.x));
            double delta = 0.0;
            if (*inputs.iter <= 1) {
                if (*inputs.udelta0 <= 0.0) {
                    delta = std::max(*inputs.delmin, 0.1 * std::max(1.0, xnorm_full));
                } else {
                    delta = *inputs.udelta0;
                }
            } else {
                delta = std::max(*inputs.delmin, 10.0 * std::sqrt(std::max(sts, 0.0)));
            }

            double acgeps = 0.0;
            double bcgeps = 0.0;
            if (inputs.cg_schedule_seed_valid) {
                acgeps = inputs.cg_schedule_acgeps;
                bcgeps = inputs.cg_schedule_bcgeps;
            } else {
                gp_ieee_signal1_cpp(
                    *gpsupn_after,
                    &acgeps,
                    &bcgeps,
                    *inputs.cgepsf,
                    *inputs.cgepsi,
                    *inputs.cggpnf
                );
            }

            const double epsgpen2 = (*inputs.epsgpen) * (*inputs.epsgpen);
            const double gpeucn20 =
                inputs.cg_schedule_seed_valid ? inputs.cg_schedule_gpeucn20 : *gpeucn2_after;
            const double gpsupn0 =
                inputs.cg_schedule_seed_valid ? inputs.cg_schedule_gpsupn0 : *gpsupn_after;
            double kappa = 0.0;
            double cgeps = *inputs.cgepsf;
            int cgmaxit = 0;
            gp_ieee_signal2_cpp(
                &cgmaxit,
                nind_work,
                *inputs.nearlyq,
                *inputs.ucgmaxit,
                *inputs.cgscre,
                &kappa,
                *gpeucn2_after,
                gpeucn20,
                epsgpen2,
                *inputs.epsgpsn,
                &cgeps,
                acgeps,
                bcgeps,
                *inputs.cgepsf,
                *inputs.cgepsi,
                *gpsupn_after,
                gpsupn0
            );

            std::vector<double> cg_s(n_val, 0.0);
            std::vector<double> cg_w(n_val, 0.0);
            std::vector<double> cg_y(n_val, 0.0);
            std::vector<double> cg_r(n_val, 0.0);
            std::vector<double> cg_d(n_val, 0.0);
            std::vector<double> cg_sprev(n_val, 0.0);
            int cg_iter = 0;
            int rbdtype = 0;
            int rbdind = 0;
            int cg_inform = 0;
            packmol_gencan_cg_bridge(
                &nind_work, inputs.ind, inputs.n, x_work.data(), inputs.m, inputs.lambda,
                inputs.rho, g_work.data(), &delta, l_work.data(), u_work.data(), &cgeps,
                inputs.epsnqmp, inputs.maxitnqmp, &cgmaxit, inputs.nearlyq, inputs.gtype,
                inputs.htvtype, inputs.trtype, inputs.iprint, inputs.ncomp, cg_s.data(),
                &cg_iter, &rbdtype, &rbdind, &cg_inform, cg_w.data(), cg_y.data(),
                cg_r.data(), cg_d.data(), cg_sprev.data(), inputs.theta, inputs.sterel,
                inputs.steabs, inputs.epsrel, inputs.epsabs, inputs.infrel, inputs.infabs
            );
            *inputs.cgcnt += cg_iter;
            if (cg_inform < 0) {
                *inputs.inform = cg_inform;
                return;
            }

            double amax = *inputs.infabs;
            if (cg_inform == 2) {
                amax = 1.0;
            } else {
                for (int i = 0; i < nind_work; ++i) {
                    if (cg_s[i] > 0.0) {
                        const double amaxx = (u_work[i] - x_work[i]) / cg_s[i];
                        if (amaxx < amax) {
                            amax = amaxx;
                            rbdind = i + 1;
                            rbdtype = 2;
                        }
                    } else if (cg_s[i] < 0.0) {
                        const double amaxx = (l_work[i] - x_work[i]) / cg_s[i];
                        if (amaxx < amax) {
                            amax = amaxx;
                            rbdind = i + 1;
                            rbdtype = 1;
                        }
                    }
                }
            }

            int tnls_inform = 0;
            const int tnint_prev = *inputs.tnintcnt;
            const int tnexg_prev = *inputs.tnexgcnt;
            const int tnexb_prev = *inputs.tnexbcnt;
            const int fcnt_prev_tn = *inputs.fcnt;
            double f_work = *inputs.f;
            std::vector<double> xplus(n_val, 0.0);
            std::vector<double> xtmp(n_val, 0.0);
            std::vector<double> xbext(n_val, 0.0);
            (void)tnls_cpp_subset(
                &nind_work, inputs.ind, inputs.n, x_work.data(), inputs.m, inputs.lambda, inputs.rho,
                l_work.data(), u_work.data(), &f_work, g_work.data(), cg_s.data(), &amax, &rbdtype,
                &rbdind, inputs.nint, inputs.next, inputs.mininterp, inputs.maxextrap, inputs.fmin,
                inputs.maxfc, inputs.gtype, inputs.fcnt, inputs.gcnt, inputs.tnintcnt, inputs.tnexgcnt,
                inputs.tnexbcnt, &tnls_inform, xplus.data(), xtmp.data(), xbext.data(), inputs.gamma,
                inputs.beta, inputs.sigma1, inputs.sigma2, inputs.sterel, inputs.steabs, inputs.epsrel,
                inputs.epsabs
            );
            if (tnls_inform < 0) {
                *inputs.inform = tnls_inform;
                return;
            }

            *inputs.tniter += 1;
            if (*inputs.tnintcnt > tnint_prev) {
                *inputs.tnintfe += (*inputs.fcnt - fcnt_prev_tn);
            } else if (*inputs.tnexgcnt > tnexg_prev) {
                *inputs.tnexgfe += (*inputs.fcnt - fcnt_prev_tn);
            } else if (*inputs.tnexbcnt > tnexb_prev) {
                *inputs.tnexbfe += (*inputs.fcnt - fcnt_prev_tn);
            } else {
                *inputs.tnstpcnt += 1;
            }
            *inputs.tnfcnt += (*inputs.fcnt - fcnt_prev_tn);

            expand_inplace(nind_work, inputs.ind, x_work.data());
            expand_inplace(nind_work, inputs.ind, g_work.data());
            expand_inplace(nind_work, inputs.ind, l_work.data());
            expand_inplace(nind_work, inputs.ind, u_work.data());

            for (int i = 0; i < n_val; ++i) {
                if (x_work[i] <= l_work[i] + std::max((*inputs.epsrel) * std::abs(l_work[i]), *inputs.epsabs)) {
                    x_work[i] = l_work[i];
                } else if (x_work[i] >= u_work[i] - std::max((*inputs.epsrel) * std::abs(u_work[i]), *inputs.epsabs)) {
                    x_work[i] = u_work[i];
                }
            }

                line_inform = tnls_inform;
                if (tnls_inform == 6) {
                    *inputs.spgiter += 1;
                    double lamspg =
                        std::max(1.0, xnorm_full) / std::sqrt(std::max(*gpeucn2_after, 1.0e-30));
                    lamspg = std::min(*inputs.lspgma, std::max(*inputs.lspgmi, lamspg));
                const int fcnt_prev_spg = *inputs.fcnt;
                int spg_line_inform = 0;
                spgls_cpp(
                    inputs.n, x_work.data(), inputs.m, inputs.lambda, inputs.rho, &f_work,
                    g_work.data(), l_work.data(), u_work.data(), &lamspg, inputs.nint,
                    inputs.mininterp, inputs.fmin, inputs.maxfc, inputs.fcnt, &spg_line_inform,
                    xtrial_work.data(), d_work.data(), inputs.gamma, inputs.sigma1, inputs.sigma2,
                    inputs.epsrel, inputs.epsabs
                );
                *inputs.spgfcnt += (*inputs.fcnt - fcnt_prev_spg);
                if (spg_line_inform < 0) {
                    *inputs.inform = spg_line_inform;
                    return;
                }
                int grad_inform = 0;
                eval_gradient_full_cpp(
                    inputs.n, x_work.data(), inputs.m, inputs.lambda, inputs.rho, inputs.gtype,
                    g_work.data(), inputs.sterel, inputs.steabs, &grad_inform
                );
                *inputs.gcnt += 1;
                if (grad_inform < 0) {
                    *inputs.inform = grad_inform;
                    return;
                }
                line_inform = spg_line_inform;
            }

            for (int i = 0; i < n_val; ++i) {
                inputs.x[i] = x_work[i];
                inputs.g[i] = g_work[i];
            }
            *inputs.f = f_work;
        } else {
            double lamspg = 0.0;
            if (*inputs.iter <= 1 || sty <= 0.0) {
                const double xnorm_for_spg = std::sqrt(norm2_kernel(n_val, inputs.x));
                lamspg = std::max(1.0, xnorm_for_spg) / std::sqrt(std::max(*gpeucn2_after, 1.0e-30));
            } else {
                lamspg = sts / sty;
            }
            lamspg = std::min(*inputs.lspgma, std::max(*inputs.lspgmi, lamspg));
            const int fcnt_before_tail = *inputs.fcnt;
            spgls_cpp(
                inputs.n, inputs.x, inputs.m, inputs.lambda, inputs.rho, inputs.f, inputs.g,
                inputs.l, inputs.u, &lamspg, inputs.nint, inputs.mininterp, inputs.fmin,
                inputs.maxfc, inputs.fcnt, &line_inform, xtrial_work.data(), d_work.data(),
                inputs.gamma, inputs.sigma1, inputs.sigma2, inputs.epsrel, inputs.epsabs
            );
            *inputs.spgfcnt += (*inputs.fcnt - fcnt_before_tail);
            *inputs.spgiter += 1;

            if (line_inform < 0) {
                *inputs.inform = line_inform;
                return;
            }

            int grad_inform = 0;
            eval_gradient_full_cpp(
                inputs.n, inputs.x, inputs.m, inputs.lambda, inputs.rho, inputs.gtype, inputs.g,
                inputs.sterel, inputs.steabs, &grad_inform
            );
            *inputs.gcnt += 1;
            if (grad_inform < 0) {
                *inputs.inform = grad_inform;
                return;
            }
        }

        project_state_and_metrics_cpp(
            n_val, inputs.x, inputs.g, inputs.l, inputs.u, *inputs.epsrel, *inputs.epsabs,
            inputs.ind, nind_after, gpsupn_after, gpeucn2_after
        );
        sts = 0.0;
        sty = 0.0;
        if (inputs.s_vec != nullptr && inputs.y_vec != nullptr) {
            for (int i = 0; i < n_val; ++i) {
                inputs.s_vec[i] = inputs.x[i] - x_prev[i];
                inputs.y_vec[i] = inputs.g[i] - g_prev[i];
                sts += inputs.s_vec[i] * inputs.s_vec[i];
                sty += inputs.s_vec[i] * inputs.y_vec[i];
            }
        }
        *inputs.gpeucn2 = *gpeucn2_after;
        *inputs.gpsupn = *gpsupn_after;
        if (*inputs.maxitngp > 0) {
            inputs.lastgpns[*inputs.iter % (*inputs.maxitngp)] = *gpeucn2_after;
        }

        int precision_after_retry_flag = 0;
        *precision_after = packmolprecision_cpp(
            inputs.n, inputs.x, inputs.m, inputs.lambda, inputs.rho, &precision_after_retry_flag
        );
        if (precision_after_retry_flag < 0) {
            *inputs.inform = precision_after_retry_flag;
            return;
        }
        *post_inform = evaluate_post_step_inform_cpp(
            *precision_after,
            line_inform,
            *gpeucn2_after,
            *inputs.epsgpen,
            *gpsupn_after,
            *inputs.epsgpsn,
            f_before_retry,
            *inputs.f,
            *inputs.epsnfp,
            *inputs.maxitnfp,
            *inputs.infabs,
            inputs.lastgpns,
            *inputs.maxitngp,
            *inputs.fmin,
            *inputs.iter,
            *inputs.maxit,
            *inputs.fcnt,
            *inputs.maxfc,
            progress_fprev,
            progress_bestprog,
            progress_itnfp
        );
        if (*post_inform >= 0) {
            break;
        }
        if (*inputs.fcnt >= *inputs.maxfc || *inputs.iter >= *inputs.maxit) {
            break;
        }
    }
}

struct ContinuationProgressState {
    double fprev;
    double bestprog;
    int itnfp;
};

ContinuationProgressState make_continuation_progress_state_cpp(
    const double f_value,
    const double infabs_value,
    const double* progress_fprev_seed,
    const double* progress_bestprog_seed,
    const int* progress_itnfp_seed
) {
    ContinuationProgressState state{
        f_value,
        std::max(infabs_value - f_value, 0.0),
        0
    };
    if (progress_fprev_seed != nullptr) {
        state.fprev = *progress_fprev_seed;
    }
    if (progress_bestprog_seed != nullptr) {
        state.bestprog = *progress_bestprog_seed;
    }
    if (progress_itnfp_seed != nullptr) {
        state.itnfp = *progress_itnfp_seed;
    }
    return state;
}

void emit_blocked_tail_debug_cpp(
    const char* fallback_reason,
    const int post_inform,
    const bool precision_after,
    const int nind_after,
    const int iter_value,
    const int fcnt_value,
    const double gpsupn_after,
    const double gpeucn2_after,
    const int continuation_budget_used,
    const int continuation_steps_done,
    const int maxfc_value
) {
    if (!gencan_debug_enabled()) {
        return;
    }
    std::fprintf(
        stderr,
        "[gencan-cpp-tail-blocked] reason=%s mode=%d post_inform=%d precision=%d nind=%d iter=%d fcnt=%d gpsupn=%.16e gpeucn2=%.16e\n",
        fallback_reason,
        static_cast<int>(active_impl_mode()),
        post_inform,
        precision_after ? 1 : 0,
        nind_after,
        iter_value,
        fcnt_value,
        gpsupn_after,
        gpeucn2_after
    );
    std::fprintf(
        stderr,
        "[gencan-cpp-tail-blocked-detail] reason=%s mode=%d budget=%d steps=%d maxfc=%d fcnt=%d\n",
        fallback_reason,
        static_cast<int>(active_impl_mode()),
        continuation_budget_used,
        continuation_steps_done,
        maxfc_value,
        fcnt_value
    );
}

void execute_tn_post_blocked_tail_cpp(
    const char* fallback_reason,
    const int n_val,
    const TnPostContinuationInputs& inputs,
    const double* progress_fprev_seed,
    const double* progress_bestprog_seed,
    const int* progress_itnfp_seed
) {
    double gpsupn_after = 0.0;
    double gpeucn2_after = 0.0;
    int nind_after = 0;
    project_state_and_metrics_cpp(
        n_val, inputs.x, inputs.g, inputs.l, inputs.u, *inputs.epsrel, *inputs.epsabs,
        inputs.ind, &nind_after, &gpsupn_after, &gpeucn2_after
    );
    *inputs.gpeucn2 = gpeucn2_after;
    *inputs.gpsupn = gpsupn_after;
    if (*inputs.maxitngp > 0) {
        inputs.lastgpns[*inputs.iter % (*inputs.maxitngp)] = gpeucn2_after;
    }

    int post_inform = -1;
    bool precision_after = false;
    int continuation_budget_used = 0;
    int continuation_steps_done = 0;
    if (*inputs.fcnt < *inputs.maxfc) {
        if (*inputs.iter <= 0) {
            *inputs.iter = 1;
        }
        ContinuationProgressState progress_state = make_continuation_progress_state_cpp(
            *inputs.f, *inputs.infabs, progress_fprev_seed, progress_bestprog_seed,
            progress_itnfp_seed
        );
        continue_tn_post_blocked_tail_cpp(
            n_val,
            inputs,
            &precision_after,
            &gpsupn_after,
            &gpeucn2_after,
            &nind_after,
            &post_inform,
            &continuation_budget_used,
            &continuation_steps_done,
            &progress_state.fprev,
            &progress_state.bestprog,
            &progress_state.itnfp
        );
    } else {
        int precision_after_flag = 0;
        precision_after = packmolprecision_cpp(
            inputs.n, inputs.x, inputs.m, inputs.lambda, inputs.rho, &precision_after_flag
        );
        if (precision_after_flag < 0) {
            *inputs.inform = precision_after_flag;
            return;
        }
    }

    emit_blocked_tail_debug_cpp(
        fallback_reason,
        post_inform,
        precision_after,
        nind_after,
        *inputs.iter,
        *inputs.fcnt,
        gpsupn_after,
        gpeucn2_after,
        continuation_budget_used,
        continuation_steps_done,
        *inputs.maxfc
    );
    *inputs.inform = finalize_blocked_tail_inform_cpp(
        post_inform,
        precision_after,
        gpeucn2_after,
        *inputs.epsgpen,
        gpsupn_after,
        *inputs.epsgpsn,
        *inputs.iter,
        *inputs.maxit,
        *inputs.fcnt,
        *inputs.maxfc
    );
}

void execute_blocked_tail_cpp(
    const char* fallback_reason,
    const int n_val,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const int* gtype,
    double* g,
    const double* l,
    const double* u,
    const double* epsrel,
    const double* epsabs,
    int* ind,
    double* gpeucn2,
    double* gpsupn,
    const int* maxitngp,
    double* lastgpns,
    int* iter,
    const int* maxfc,
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
    double* s_vec,
    double* y_vec,
    double* w_vec,
    const double* udelta0,
    const int* ucgmaxit,
    const int* cgscre,
    const double* cggpnf,
    const double* cgepsi,
    const double* cgepsf,
    const double* epsnqmp,
    const int* maxitnqmp,
    const bool* nearlyq,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    const double* eta,
    const double* delmin,
    const double* lspgma,
    const double* lspgmi,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    double* f,
    const double* theta,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsgpen,
    const double* epsgpsn,
    const int* maxitnfp,
    const double* epsnfp,
    const double* infabs,
    const int* maxit,
    int* inform,
    const double* progress_fprev_seed,
    const double* progress_bestprog_seed,
    const int* progress_itnfp_seed
) {
    const std::string blocked_reason(fallback_reason);
    (void)w_vec;
    (void)ncomp;
    int grad_inform = 0;
    eval_gradient_full_cpp(n, x, m, lambda, rho, gtype, g, sterel, steabs, &grad_inform);
    *gcnt += 1;
    if (grad_inform < 0) {
        *inform = grad_inform;
        return;
    }

    double gpsupn_after = 0.0;
    double gpeucn2_after = 0.0;
    int nind_after = 0;
    project_state_and_metrics_cpp(
        n_val, x, g, l, u, *epsrel, *epsabs, ind, &nind_after, &gpsupn_after, &gpeucn2_after
    );
    *gpeucn2 = gpeucn2_after;
    *gpsupn = gpsupn_after;
    if (*maxitngp > 0) {
        lastgpns[*iter % (*maxitngp)] = gpeucn2_after;
    }

    const bool can_continue_from_post =
        blocked_reason == "spg_post_nonterminal" ||
        blocked_reason == "tn_no_free_variables";
    int precision_after_flag = 0;
    bool precision_after = packmolprecision_cpp(
        n, x, m, lambda, rho, &precision_after_flag
    );
    if (precision_after_flag < 0) {
        *inform = precision_after_flag;
        return;
    }
    int post_inform = -1;
    int continuation_budget_used = 0;
    int continuation_steps_done = 0;
    if (can_continue_from_post && *fcnt < *maxfc) {
        if (*iter <= 0) {
            *iter = 1;
        }
        ContinuationProgressState progress_state = make_continuation_progress_state_cpp(
            *f, *infabs, progress_fprev_seed, progress_bestprog_seed, progress_itnfp_seed
        );
        post_inform = evaluate_post_step_inform_cpp(
            precision_after,
            0,
            gpeucn2_after,
            *epsgpen,
            gpsupn_after,
            *epsgpsn,
            *f,
            *f,
            *epsnfp,
            *maxitnfp,
            *infabs,
            lastgpns,
            *maxitngp,
            *fmin,
            *iter,
            *maxit,
            *fcnt,
            *maxfc,
            &progress_state.fprev,
            &progress_state.bestprog,
            &progress_state.itnfp
        );
        if (post_inform >= 0) {
            continuation_budget_used = 0;
            continuation_steps_done = 0;
        } else {
            const int continuation_budget = block_tail_continuation_steps();
            const bool unlimited_continuation = continuation_budget <= 0;
            continuation_budget_used = continuation_budget;
            std::vector<double> xtrial_work(n_val, 0.0);
            std::vector<double> d_work(n_val, 0.0);
            std::vector<double> x_prev(n_val, 0.0);
            std::vector<double> g_prev(n_val, 0.0);
            double sts = 0.0;
            double sty = 0.0;
            if (s_vec != nullptr && y_vec != nullptr) {
                for (int i = 0; i < n_val; ++i) {
                    sts += s_vec[i] * s_vec[i];
                    sty += y_vec[i] * y_vec[i];
                }
            }
            while (unlimited_continuation || continuation_steps_done < continuation_budget) {
                continuation_steps_done += 1;
                *iter += 1;
                const double f_before_retry = *f;
                for (int i = 0; i < n_val; ++i) {
                    x_prev[i] = x[i];
                    g_prev[i] = g[i];
                }
                int line_inform = 0;
                double gieucn2_after = 0.0;
                for (int i = 0; i < n_val; ++i) {
                    const double xpg = x[i] - g[i];
                    const double gpi = std::min(u[i], std::max(l[i], xpg)) - x[i];
                    if (x[i] > l[i] && x[i] < u[i]) {
                        gieucn2_after += gpi * gpi;
                    }
                }
                const double ometa2 = (1.0 - *eta) * (1.0 - *eta);
                const bool use_spg_step =
                    (gpeucn2_after > 0.0 && gieucn2_after <= ometa2 * gpeucn2_after) ||
                    nind_after <= 0;

            if (!use_spg_step) {
                std::vector<double> x_work(x, x + n_val);
                std::vector<double> g_work(g, g + n_val);
                std::vector<double> l_work(l, l + n_val);
                std::vector<double> u_work(u, u + n_val);
                const int nind_work = nind_after;

                shrink_inplace(nind_work, ind, x_work.data());
                shrink_inplace(nind_work, ind, g_work.data());
                shrink_inplace(nind_work, ind, l_work.data());
                shrink_inplace(nind_work, ind, u_work.data());

                const double xnorm_full = std::sqrt(norm2_kernel(n_val, x));
                double delta = 0.0;
                if (*udelta0 <= 0.0) {
                    delta = std::max(*delmin, 0.1 * std::max(1.0, xnorm_full));
                } else {
                    delta = *udelta0;
                }

                double acgeps = 0.0;
                double bcgeps = 0.0;
                gp_ieee_signal1_cpp(gpsupn_after, &acgeps, &bcgeps, *cgepsf, *cgepsi, *cggpnf);

                const double epsgpen2 = (*epsgpen) * (*epsgpen);
                const double gpeucn20 = gpeucn2_after;
                const double gpsupn0 = gpsupn_after;
                double kappa = 0.0;
                double cgeps = *cgepsf;
                int cgmaxit = 0;
                gp_ieee_signal2_cpp(
                    &cgmaxit, nind_work, *nearlyq, *ucgmaxit, *cgscre, &kappa, gpeucn2_after,
                    gpeucn20, epsgpen2, *epsgpsn, &cgeps, acgeps, bcgeps, *cgepsf, *cgepsi,
                    gpsupn_after, gpsupn0
                );

                std::vector<double> cg_s(n_val, 0.0);
                std::vector<double> cg_w(n_val, 0.0);
                std::vector<double> cg_y(n_val, 0.0);
                std::vector<double> cg_r(n_val, 0.0);
                std::vector<double> cg_d(n_val, 0.0);
                std::vector<double> cg_sprev(n_val, 0.0);
                int cg_iter = 0;
                int rbdtype = 0;
                int rbdind = 0;
                int cg_inform = 0;
                cg_cpp_full(
                    &nind_work, ind, n, x_work.data(), m, lambda, rho, g_work.data(), &delta,
                    l_work.data(), u_work.data(), &cgeps, epsnqmp, maxitnqmp, &cgmaxit, nearlyq,
                    gtype, htvtype, trtype, iprint, theta, sterel, steabs, epsrel, epsabs, infabs,
                    cg_s.data(), cg_w.data(), cg_y.data(), cg_r.data(), cg_d.data(), cg_sprev.data(),
                    &cg_iter, &rbdtype, &rbdind, &cg_inform
                );
                *cgcnt += cg_iter;
                if (cg_inform < 0) {
                    *inform = cg_inform;
                    return;
                }

                double amax = *infabs;
                if (cg_inform == 2) {
                    amax = 1.0;
                } else {
                    for (int i = 0; i < nind_work; ++i) {
                        if (cg_s[i] > 0.0) {
                            const double amaxx = (u_work[i] - x_work[i]) / cg_s[i];
                            if (amaxx < amax) {
                                amax = amaxx;
                                rbdind = i + 1;
                                rbdtype = 2;
                            }
                        } else if (cg_s[i] < 0.0) {
                            const double amaxx = (l_work[i] - x_work[i]) / cg_s[i];
                            if (amaxx < amax) {
                                amax = amaxx;
                                rbdind = i + 1;
                                rbdtype = 1;
                            }
                        }
                    }
                }

                int tnls_inform = 0;
                const int tnint_prev = *tnintcnt;
                const int tnexg_prev = *tnexgcnt;
                const int tnexb_prev = *tnexbcnt;
                const int fcnt_prev_tn = *fcnt;
                double f_work = *f;
                std::vector<double> xplus(n_val, 0.0);
                std::vector<double> xtmp(n_val, 0.0);
                std::vector<double> xbext(n_val, 0.0);
                (void)tnls_cpp_subset(
                    &nind_work, ind, n, x_work.data(), m, lambda, rho, l_work.data(), u_work.data(),
                    &f_work, g_work.data(), cg_s.data(), &amax, &rbdtype, &rbdind, nint, next, mininterp,
                    maxextrap, fmin, maxfc, gtype, fcnt, gcnt, tnintcnt, tnexgcnt, tnexbcnt, &tnls_inform,
                    xplus.data(), xtmp.data(), xbext.data(), gamma, beta, sigma1, sigma2, sterel, steabs,
                    epsrel, epsabs
                );
                if (tnls_inform < 0) {
                    *inform = tnls_inform;
                    return;
                }

                *tniter += 1;
                if (*tnintcnt > tnint_prev) {
                    *tnintfe += (*fcnt - fcnt_prev_tn);
                } else if (*tnexgcnt > tnexg_prev) {
                    *tnexgfe += (*fcnt - fcnt_prev_tn);
                } else if (*tnexbcnt > tnexb_prev) {
                    *tnexbfe += (*fcnt - fcnt_prev_tn);
                } else {
                    *tnstpcnt += 1;
                }
                *tnfcnt += (*fcnt - fcnt_prev_tn);

                expand_inplace(nind_work, ind, x_work.data());
                expand_inplace(nind_work, ind, g_work.data());
                expand_inplace(nind_work, ind, l_work.data());
                expand_inplace(nind_work, ind, u_work.data());

                for (int i = 0; i < n_val; ++i) {
                    if (x_work[i] <= l_work[i] + std::max((*epsrel) * std::abs(l_work[i]), *epsabs)) {
                        x_work[i] = l_work[i];
                    } else if (x_work[i] >= u_work[i] - std::max((*epsrel) * std::abs(u_work[i]), *epsabs)) {
                        x_work[i] = u_work[i];
                    }
                }

                line_inform = tnls_inform;
                if (tnls_inform == 6) {
                    *spgiter += 1;
                    double lamspg = std::max(1.0, xnorm_full) / std::sqrt(std::max(gpeucn2_after, 1.0e-30));
                    lamspg = std::min(*lspgma, std::max(*lspgmi, lamspg));
                    const int fcnt_prev_spg = *fcnt;
                    int spg_line_inform = 0;
                    spgls_cpp(
                        n, x_work.data(), m, lambda, rho, &f_work, g_work.data(), l_work.data(), u_work.data(),
                        &lamspg, nint, mininterp, fmin, maxfc, fcnt, &spg_line_inform, xtrial_work.data(),
                        d_work.data(), gamma, sigma1, sigma2, epsrel, epsabs
                    );
                    *spgfcnt += (*fcnt - fcnt_prev_spg);
                    if (spg_line_inform < 0) {
                        *inform = spg_line_inform;
                        return;
                    }
                    grad_inform = 0;
                    eval_gradient_full_cpp(
                        n, x_work.data(), m, lambda, rho, gtype, g_work.data(), sterel, steabs, &grad_inform
                    );
                    *gcnt += 1;
                    if (grad_inform < 0) {
                        *inform = grad_inform;
                        return;
                    }
                    line_inform = spg_line_inform;
                }

                for (int i = 0; i < n_val; ++i) {
                    x[i] = x_work[i];
                    g[i] = g_work[i];
                }
                *f = f_work;

                grad_inform = 0;
                eval_gradient_full_cpp(n, x, m, lambda, rho, gtype, g, sterel, steabs, &grad_inform);
                *gcnt += 1;
                if (grad_inform < 0) {
                    *inform = grad_inform;
                    return;
                }
            } else {
                double lamspg = 0.0;
                if (*iter <= 1 || sty <= 0.0) {
                    const double xnorm_for_spg = std::sqrt(norm2_kernel(n_val, x));
                    lamspg = std::max(1.0, xnorm_for_spg) / std::sqrt(std::max(gpeucn2_after, 1.0e-30));
                } else {
                    lamspg = sts / sty;
                }
                lamspg = std::min(*lspgma, std::max(*lspgmi, lamspg));
                const int fcnt_before_tail = *fcnt;
                spgls_cpp(
                    n, x, m, lambda, rho, f, g, l, u, &lamspg, nint, mininterp, fmin, maxfc, fcnt,
                    &line_inform, xtrial_work.data(), d_work.data(), gamma, sigma1, sigma2, epsrel, epsabs
                );
                *spgfcnt += (*fcnt - fcnt_before_tail);
                *spgiter += 1;

                if (line_inform < 0) {
                    *inform = line_inform;
                    return;
                }

                grad_inform = 0;
                eval_gradient_full_cpp(n, x, m, lambda, rho, gtype, g, sterel, steabs, &grad_inform);
                *gcnt += 1;
                if (grad_inform < 0) {
                    *inform = grad_inform;
                    return;
                }
            }

            project_state_and_metrics_cpp(
                n_val, x, g, l, u, *epsrel, *epsabs, ind, &nind_after, &gpsupn_after, &gpeucn2_after
            );
            sts = 0.0;
            sty = 0.0;
            if (s_vec != nullptr && y_vec != nullptr) {
                for (int i = 0; i < n_val; ++i) {
                    s_vec[i] = x[i] - x_prev[i];
                    y_vec[i] = g[i] - g_prev[i];
                    sts += s_vec[i] * s_vec[i];
                    sty += s_vec[i] * y_vec[i];
                }
            }
            *gpeucn2 = gpeucn2_after;
            *gpsupn = gpsupn_after;
            if (*maxitngp > 0) {
                lastgpns[*iter % (*maxitngp)] = gpeucn2_after;
            }

            int precision_after_retry_flag = 0;
            precision_after = packmolprecision_cpp(
                n, x, m, lambda, rho, &precision_after_retry_flag
            );
            if (precision_after_retry_flag < 0) {
                *inform = precision_after_retry_flag;
                return;
            }
            post_inform = evaluate_post_step_inform_cpp(
                precision_after,
                line_inform,
                gpeucn2_after,
                *epsgpen,
                gpsupn_after,
                *epsgpsn,
                f_before_retry,
                *f,
                *epsnfp,
                *maxitnfp,
                *infabs,
                lastgpns,
                *maxitngp,
                *fmin,
                *iter,
                *maxit,
                *fcnt,
                *maxfc,
                &progress_state.fprev,
                &progress_state.bestprog,
                &progress_state.itnfp
            );
            if (post_inform >= 0) {
                break;
            }
            if (*fcnt >= *maxfc || *iter >= *maxit) {
                break;
            }
        }
    }
    }

    emit_blocked_tail_debug_cpp(
        fallback_reason,
        post_inform,
        precision_after,
        nind_after,
        *iter,
        *fcnt,
        gpsupn_after,
        gpeucn2_after,
        continuation_budget_used,
        continuation_steps_done,
        *maxfc
    );
    *inform = finalize_blocked_tail_inform_cpp(
        post_inform,
        precision_after,
        gpeucn2_after,
        *epsgpen,
        gpsupn_after,
        *epsgpsn,
        *iter,
        *maxit,
        *fcnt,
        *maxfc
    );
}

extern "C" void packmol_gencan_fortran_c(
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

extern "C" void packmol_gencan_gencan_bridge(
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

bool try_spg_post_cpp_replay_tail_cpp(
    const char* fallback_reason,
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
) {
    if (fallback_reason == nullptr || std::string(fallback_reason) != "spg_post_nonterminal") {
        return false;
    }
    if (active_impl_mode() != GencanImplMode::kCpp) {
        return false;
    }
    const bool ab_compare_spg_post_replay = ab_compare_enabled();
    if (!ab_compare_spg_post_replay &&
        (!handoff_cpp_tail_enabled() || !handoff_cpp_tail_allow_spg_post_enabled())) {
        return false;
    }
    if (!ab_compare_spg_post_replay &&
        init1_phase_active_cpp() &&
        !handoff_cpp_tail_allow_init1_enabled()) {
        return false;
    }
    const int max_depth = spg_post_cpp_replay_max_depth();
    if (g_spg_post_cpp_replay_depth >= max_depth) {
        return false;
    }
    g_spg_post_cpp_replay_depth += 1;
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp-spg-post-cpp-replay] depth=%d max_depth=%d mode=%d\n",
            g_spg_post_cpp_replay_depth,
            max_depth,
            static_cast<int>(active_impl_mode())
        );
    }
    packmol_gencan_gencan_bridge(
        n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit,
        maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq,
        nint, next, mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2,
        gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt,
        tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform, s, y, d, ind, lastgpns, w, eta,
        delmin, lspgma, lspgmi, theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs,
        infrel, infabs
    );
    if (*inform >= 0) {
        apply_spg_post_cpp_replay_bookkeeping(iter, fcnt, gcnt);
        if (gencan_debug_enabled()) {
            std::fprintf(
                stderr,
                "[gencan-cpp-spg-post-cpp-replay-bookkeeping] iter=%d fcnt=%d gcnt=%d\n",
                *iter,
                *fcnt,
                *gcnt
            );
        }
    }
    g_spg_post_cpp_replay_depth -= 1;
    return true;
}

bool try_tn_post_cpp_replay_tail_cpp(
    const char* fallback_reason,
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
) {
    const std::string reason = fallback_reason == nullptr ? "" : std::string(fallback_reason);
    if (fallback_reason == nullptr || !is_tn_tail_reason_cpp(reason)) {
        return false;
    }
    const bool nested_spg_replay = g_spg_post_cpp_replay_depth > 0;
    if (nested_spg_replay && !handoff_cpp_tail_allow_tn_post_explicitly_enabled()) {
        return false;
    }
    const bool replay_enabled = tn_post_cpp_replay_tail_enabled();
    const bool cpp_mode = active_impl_mode() == GencanImplMode::kCpp && !ab_compare_enabled();
    const int max_depth = tn_post_cpp_replay_max_depth();
    const bool depth_ok = g_tn_post_cpp_replay_depth < max_depth;
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp-tn-post-cpp-replay-gate] reason=%s mode=%d replay_enabled=%d cpp_mode=%d depth=%d max_depth=%d\n",
            reason.c_str(),
            static_cast<int>(active_impl_mode()),
            replay_enabled ? 1 : 0,
            cpp_mode ? 1 : 0,
            g_tn_post_cpp_replay_depth,
            max_depth
        );
    }
    if (!replay_enabled || !cpp_mode || !depth_ok) {
        return false;
    }
    g_tn_post_cpp_replay_depth += 1;
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp-tn-post-cpp-replay] depth=%d max_depth=%d mode=%d\n",
            g_tn_post_cpp_replay_depth,
            max_depth,
            static_cast<int>(active_impl_mode())
        );
    }
    packmol_gencan_gencan_bridge(
        n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit,
        maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq,
        nint, next, mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2,
        gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt,
        tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform, s, y, d, ind, lastgpns, w, eta,
        delmin, lspgma, lspgmi, theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs,
        infrel, infabs
    );
    g_tn_post_cpp_replay_depth -= 1;
    return true;
}

void emit_tn_shadow_report_cpp(
    const int n_val,
    const double* x,
    const double* g,
    const std::vector<double>& x_before,
    const std::vector<double>& g_before,
    const double f_before,
    const int inform_before,
    const int iter_before,
    const int fcnt_before,
    const int gcnt_before,
    const int cgcnt_before,
    const int spgiter_before,
    const int tniter_before,
    const double f_after,
    const int inform_after,
    const int iter_after,
    const int fcnt_after,
    const int gcnt_after,
    const int cgcnt_after,
    const int spgiter_after,
    const int tniter_after
) {
    double max_abs_dx = 0.0;
    double max_abs_dg = 0.0;
    for (int i = 0; i < n_val; ++i) {
        max_abs_dx = std::max(max_abs_dx, std::abs(x[i] - x_before[i]));
        max_abs_dg = std::max(max_abs_dg, std::abs(g[i] - g_before[i]));
    }
    std::fprintf(
        stderr,
        "[gencan-tn-post-shadow] mode=%d inform:%d->%d iter:+%d fcnt:+%d gcnt:+%d cgcnt:+%d spgiter:+%d tniter:+%d df=%.16e max|dx|=%.16e max|dg|=%.16e\n",
        static_cast<int>(active_impl_mode()),
        inform_before,
        inform_after,
        iter_after - iter_before,
        fcnt_after - fcnt_before,
        gcnt_after - gcnt_before,
        cgcnt_after - cgcnt_before,
        spgiter_after - spgiter_before,
        tniter_after - tniter_before,
        f_after - f_before,
        max_abs_dx,
        max_abs_dg
    );
}

struct TnShadowSnapshot {
    bool enabled = false;
    std::vector<double> x_before;
    std::vector<double> g_before;
    double f_before = 0.0;
    int inform_before = 0;
    int iter_before = 0;
    int fcnt_before = 0;
    int gcnt_before = 0;
    int cgcnt_before = 0;
    int spgiter_before = 0;
    int tniter_before = 0;
};

TnShadowSnapshot capture_tn_shadow_snapshot_cpp(
    const bool enabled,
    const int n_val,
    const double* x,
    const double* g,
    const double* f,
    const int* inform,
    const int* iter,
    const int* fcnt,
    const int* gcnt,
    const int* cgcnt,
    const int* spgiter,
    const int* tniter
) {
    TnShadowSnapshot snapshot;
    snapshot.enabled = enabled;
    if (!enabled) {
        return snapshot;
    }
    snapshot.x_before.assign(x, x + n_val);
    snapshot.g_before.assign(g, g + n_val);
    snapshot.f_before = *f;
    snapshot.inform_before = *inform;
    snapshot.iter_before = *iter;
    snapshot.fcnt_before = *fcnt;
    snapshot.gcnt_before = *gcnt;
    snapshot.cgcnt_before = *cgcnt;
    snapshot.spgiter_before = *spgiter;
    snapshot.tniter_before = *tniter;
    return snapshot;
}

void maybe_emit_tn_shadow_report_cpp(
    const TnShadowSnapshot& snapshot,
    const int n_val,
    const double* x,
    const double* g,
    const double f_after,
    const int inform_after,
    const int iter_after,
    const int fcnt_after,
    const int gcnt_after,
    const int cgcnt_after,
    const int spgiter_after,
    const int tniter_after
) {
    if (!snapshot.enabled) {
        return;
    }
    emit_tn_shadow_report_cpp(
        n_val,
        x,
        g,
        snapshot.x_before,
        snapshot.g_before,
        snapshot.f_before,
        snapshot.inform_before,
        snapshot.iter_before,
        snapshot.fcnt_before,
        snapshot.gcnt_before,
        snapshot.cgcnt_before,
        snapshot.spgiter_before,
        snapshot.tniter_before,
        f_after,
        inform_after,
        iter_after,
        fcnt_after,
        gcnt_after,
        cgcnt_after,
        spgiter_after,
        tniter_after
    );
}

void continue_spg_post_nonterminal_cpp(
    const int n_val,
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
    const double* infabs,
    const double progress_fprev,
    const double progress_bestprog,
    const int progress_itnfp
) {
    if (try_spg_post_cpp_replay_tail_cpp(
            "spg_post_nonterminal", n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp,
            maxitngp, fmin, maxit, maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf,
            epsnqmp, maxitnqmp, nearlyq, nint, next, mininterp, maxextrap, gtype, htvtype,
            trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter, fcnt, gcnt, cgcnt, spgiter,
            spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe,
            tnexbfe, inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta,
            gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
        )) {
        return;
    }
    execute_blocked_tail_cpp(
        "spg_post_nonterminal", n_val, n, x, m, lambda, rho, gtype, g, l, u, epsrel, epsabs, ind,
        gpeucn2, gpsupn, maxitngp, lastgpns, iter, maxfc, fcnt, gcnt, cgcnt, spgiter, spgfcnt,
        tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe,
        s, y, w, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq,
        htvtype, trtype, iprint, ncomp, eta, delmin, lspgma, lspgmi, nint, next, mininterp,
        maxextrap, fmin, f, theta, gamma, beta, sigma1, sigma2, sterel, steabs,
        epsgpen, epsgpsn, maxitnfp, epsnfp, infabs, maxit, inform,
        &progress_fprev, &progress_bestprog, &progress_itnfp
    );
}

void continue_tn_tail_reason_cpp(
    const char* fallback_reason,
    const int n_val,
    const TnShadowSnapshot& tn_shadow,
    const bool fallback_progress_seed_valid,
    const double fallback_progress_fprev,
    const double fallback_progress_bestprog,
    const int fallback_progress_itnfp,
    const bool fallback_cg_schedule_seed_valid,
    const double fallback_cg_schedule_acgeps,
    const double fallback_cg_schedule_bcgeps,
    const double fallback_cg_schedule_gpeucn20,
    const double fallback_cg_schedule_gpsupn0,
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
) {
    if (try_tn_post_cpp_replay_tail_cpp(
            fallback_reason, n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp,
            maxitngp, fmin, maxit, maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf,
            epsnqmp, maxitnqmp, nearlyq, nint, next, mininterp, maxextrap, gtype, htvtype,
            trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter, fcnt, gcnt, cgcnt, spgiter,
            spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe,
            tnexbfe, inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta,
            gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
        )) {
        maybe_emit_tn_shadow_report_cpp(
            tn_shadow, n_val, x, g, *f, *inform, *iter, *fcnt, *gcnt, *cgcnt, *spgiter, *tniter
        );
        return;
    }
    const TnPostContinuationInputs continuation_inputs{
        n, x, m, lambda, rho, gtype, g, l, u, epsrel, epsabs, ind, gpeucn2, gpsupn,
        maxitngp, lastgpns, iter, maxfc, fcnt, gcnt, cgcnt, spgiter, spgfcnt,
        tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe,
        tnexbfe, s, y, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp,
        maxitnqmp, nearlyq, htvtype, trtype, iprint, ncomp, eta, delmin, lspgma,
        lspgmi, nint, next, mininterp, maxextrap, fmin, f, theta, gamma, beta,
        sigma1, sigma2, sterel, steabs, epsgpen, epsgpsn, maxitnfp, epsnfp, infrel,
        infabs, maxit, fallback_cg_schedule_seed_valid, fallback_cg_schedule_acgeps,
        fallback_cg_schedule_bcgeps, fallback_cg_schedule_gpeucn20, fallback_cg_schedule_gpsupn0,
        inform
    };
    execute_tn_post_blocked_tail_cpp(
        fallback_reason,
        n_val,
        continuation_inputs,
        fallback_progress_seed_valid ? &fallback_progress_fprev : nullptr,
        fallback_progress_seed_valid ? &fallback_progress_bestprog : nullptr,
        fallback_progress_seed_valid ? &fallback_progress_itnfp : nullptr
    );
    maybe_emit_tn_shadow_report_cpp(
        tn_shadow, n_val, x, g, *f, *inform, *iter, *fcnt, *gcnt, *cgcnt, *spgiter, *tniter
    );
}

}  // namespace

extern "C" int packmol_gencan_cpp_probe(int n, const double* x, double* fx_out) {
    if (n < 0 || x == nullptr || fx_out == nullptr) {
        return -1;
    }

    // Phase B skeleton: keep behavior-neutral bridge available for mixed builds.
    *fx_out = 0.0;
    return 0;
}

extern "C" int packmol_gencan_impl_mode_c() {
    return static_cast<int>(active_impl_mode());
}

extern "C" void packmol_spgls_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    const double* g,
    const double* l,
    const double* u,
    const double* lamspg,
    const double* nint,
    const int* mininterp,
    const double* fmin,
    const int* maxfc,
    const int* iprint,
    int* fcnt,
    int* inform,
    double* xtrial,
    double* d,
    const double* gamma,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);

extern "C" void packmol_evalal_fortran_c(
    const int* n,
    const double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* flag
);

extern "C" void packmol_evalnal_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    int* flag
);

extern "C" void packmol_evalnaldiff_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* sterel,
    const double* steabs,
    int* flag
);

extern "C" void packmol_precision_state_fortran_c(
    double* precision_value,
    double* fdist_value,
    double* frest_value
);

extern "C" void packmol_computef_fortran_c(
    const int* n,
    const double* x,
    double* f
);

extern "C" void packmol_computeg_fortran_c(
    const int* n,
    const double* x,
    double* g
);

extern "C" double packmol_input_precision;
extern "C" double packmol_fdist;
extern "C" double packmol_frest;

extern "C" void packmol_evalal_fortran_c(
    const int* n,
    const double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* flag
) {
    (void)m;
    (void)lambda;
    (void)rho;
    packmol_computef_fortran_c(n, x, f);
    *flag = 0;
}

extern "C" void packmol_evalnal_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    int* flag
) {
    (void)m;
    (void)lambda;
    (void)rho;
    packmol_computeg_fortran_c(n, x, g);
    *flag = 0;
}

extern "C" void packmol_evalnaldiff_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* sterel,
    const double* steabs,
    int* flag
) {
    (void)sterel;
    (void)steabs;
    packmol_evalnal_fortran_c(n, x, m, lambda, rho, g, flag);
}

extern "C" void packmol_precision_state_fortran_c(
    double* precision_value,
    double* fdist_value,
    double* frest_value
) {
    *precision_value = packmol_input_precision;
    *fdist_value = packmol_fdist;
    *frest_value = packmol_frest;
}

extern "C" void packmol_gencan_set_init1_phase_c(
    const int* flag
) {
    g_init1_phase_active = (flag != nullptr && *flag != 0) ? 1 : 0;
}

namespace {

double dot_cpp_stable(const int n, const double* a, const double* b) {
    double sum = 0.0;
    double c = 0.0;
    for (int i = 0; i < n; ++i) {
        const double prod = a[i] * b[i];
        const double y = prod - c;
        const double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

double norm2sq_cpp_stable(const int n, const double* x) {
    double scale = 0.0;
    double ssq = 1.0;
    for (int i = 0; i < n; ++i) {
        const double ax = std::abs(x[i]);
        if (ax == 0.0) {
            continue;
        }
        if (scale < ax) {
            const double r = (scale == 0.0) ? 0.0 : (scale / ax);
            ssq = 1.0 + ssq * r * r;
            scale = ax;
        } else {
            const double r = ax / scale;
            ssq += r * r;
        }
    }
    if (scale == 0.0) {
        return 0.0;
    }
    return scale * scale * ssq;
}

double dot_cpp_legacy(const int n, const double* a, const double* b) {
    double val = 0.0;
    for (int i = 0; i < n; ++i) {
        val += a[i] * b[i];
    }
    return val;
}

double hsldnrm2_cpp_legacy(const int n, const double* x) {
    constexpr double kZero = 0.0;
    constexpr double kOne = 1.0;
    constexpr double kCutlo = 8.232e-11;
    constexpr double kCuthi = 1.304e19;
    if (n <= 0) {
        return kZero;
    }

    double sum = kZero;
    int i = 0;
    while (i < n) {
        if (std::abs(x[i]) > kCutlo) {
            break;
        }
        double xmax = kZero;
        if (x[i] == kZero) {
            i += 1;
            continue;
        }
        if (std::abs(x[i]) > kCutlo) {
            break;
        }
        xmax = std::abs(x[i]);
        while (true) {
            if (std::abs(x[i]) > kCutlo) {
                break;
            }
            if (std::abs(x[i]) > xmax) {
                sum = kOne + sum * (xmax / x[i]) * (xmax / x[i]);
                xmax = std::abs(x[i]);
            } else {
                sum += (x[i] / xmax) * (x[i] / xmax);
            }
            i += 1;
            if (i >= n) {
                return xmax * std::sqrt(sum);
            }
        }
        sum = (sum * xmax) * xmax;
        const double hitest = kCuthi / static_cast<double>(n);
        for (int j = i; j < n; ++j) {
            if (std::abs(x[j]) >= hitest) {
                i = j;
                sum = (sum / x[i]) / x[i];
                break;
            }
            sum += x[j] * x[j];
            if (j == n - 1) {
                return std::sqrt(sum);
            }
        }
    }

    const double hitest = kCuthi / static_cast<double>(n);
    for (int j = i; j < n; ++j) {
        if (std::abs(x[j]) >= hitest) {
            i = j;
            sum = (sum / x[i]) / x[i];
            double xmax = std::abs(x[i]);
            i += 1;
            while (i < n) {
                if (std::abs(x[i]) <= kCutlo) {
                    if (std::abs(x[i]) > xmax) {
                        sum = kOne + sum * (xmax / x[i]) * (xmax / x[i]);
                        xmax = std::abs(x[i]);
                    } else {
                        sum += (x[i] / xmax) * (x[i] / xmax);
                    }
                    i += 1;
                    continue;
                }
                sum = (sum * xmax) * xmax;
                break;
            }
            if (i >= n) {
                return xmax * std::sqrt(sum);
            }
            for (int k = i; k < n; ++k) {
                sum += x[k] * x[k];
            }
            return std::sqrt(sum);
        }
        sum += x[j] * x[j];
    }
    return std::sqrt(sum);
}

double dot_kernel(const int n, const double* a, const double* b) {
    if (use_cpp_numeric_kernel()) {
        return dot_cpp_stable(n, a, b);
    }
    return dot_cpp_legacy(n, a, b);
}

double norm2_kernel(const int n, const double* x) {
    if (use_cpp_numeric_kernel()) {
        return norm2sq_cpp_stable(n, x);
    }
    const double nrm = hsldnrm2_cpp_legacy(n, x);
    return nrm * nrm;
}

void vec_copy(const int n, const double* src, double* dst) {
    for (int i = 0; i < n; ++i) {
        dst[i] = src[i];
    }
}

void vec_trial_point(const int n, const double* x, const double alpha, const double* d, double* out) {
    for (int i = 0; i < n; ++i) {
        out[i] = x[i] + alpha * d[i];
    }
}

bool same_step_relative(
    const int n,
    const double alpha,
    const double* d,
    const double* ref,
    const double epsrel,
    const double epsabs
) {
    for (int i = 0; i < n; ++i) {
        if (std::abs(alpha * d[i]) > std::max(epsrel * std::abs(ref[i]), epsabs)) {
            return false;
        }
    }
    return true;
}

bool same_point_relative(
    const int n,
    const double* a,
    const double* b,
    const double* ref,
    const double epsrel,
    const double epsabs
) {
    for (int i = 0; i < n; ++i) {
        if (std::abs(a[i] - b[i]) > std::max(epsrel * std::abs(ref[i]), epsabs)) {
            return false;
        }
    }
    return true;
}

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
) {
    const int nind_val = *nind;
    const int n_val = *n;
    for (int i = nind_val; i < n_val; ++i) {
        x[i] = xc[i];
    }
    expand_inplace(nind_val, ind, x);
    packmol_evalal_fortran_c(n, x, m, lambda, rho, f, inform);
    shrink_inplace(nind_val, ind, x);
}

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
) {
    const int nind_val = *nind;
    const int n_val = *n;
    for (int i = nind_val; i < n_val; ++i) {
        x[i] = xc[i];
    }
    expand_inplace(nind_val, ind, x);
    if (*gtype == 0) {
        packmol_evalnal_fortran_c(n, x, m, lambda, rho, g, inform);
    } else {
        packmol_evalnaldiff_fortran_c(n, x, m, lambda, rho, g, sterel, steabs, inform);
    }
    shrink_inplace(nind_val, ind, x);
    shrink_inplace(nind_val, ind, g);
}

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
) {
    if (*gtype == 0) {
        packmol_evalnal_fortran_c(n, x, m, lambda, rho, g, inform);
    } else {
        packmol_evalnaldiff_fortran_c(n, x, m, lambda, rho, g, sterel, steabs, inform);
    }
}

bool packmolprecision_from_state_cpp() {
    double precision_value = 0.0;
    double fdist_value = 0.0;
    double frest_value = 0.0;
    packmol_precision_state_fortran_c(&precision_value, &fdist_value, &frest_value);
    return fdist_value < precision_value && frest_value < precision_value;
}

bool init1_phase_active_cpp() {
    return g_init1_phase_active != 0;
}

struct TailFallbackPolicy {
    bool tn_post_case;
    bool spg_post_case;
    bool tn_handoff_requested;
    bool spg_post_handoff;
    bool allow_cpp_tail_reduction;
    bool handoff_cpp_tail;
    bool handoff_cpp_tail_init1_guarded;
    bool handoff_cpp_tail_spg_guarded;
    bool suppress_fallback_marker;
    bool effective_block_tail;
};

TailFallbackPolicy make_tail_fallback_policy(const char* fallback_reason) {
    const std::string reason(fallback_reason == nullptr ? "" : fallback_reason);
    const bool tn_post_case = is_tn_tail_reason_cpp(reason);
    const bool spg_post_case = reason == "spg_post_nonterminal";
    const bool tn_handoff_requested =
        tn_post_handoff_enabled() || tn_post_handoff_safe_enabled() || tn_post_handoff_unsafe_enabled();
    const bool spg_post_handoff = spg_post_handoff_enabled() && spg_post_case;
    const bool init1_active = init1_phase_active_cpp();
    bool tn_tail_allowed = !tn_post_case || handoff_cpp_tail_allow_tn_post_enabled();
    if (tn_post_case &&
        g_spg_post_cpp_replay_depth > 0 &&
        !handoff_cpp_tail_allow_tn_post_explicitly_enabled()) {
        // The nested TN tail reached from SPG replay still drifts in state and
        // counters unless the caller explicitly opts into that experimental path.
        tn_tail_allowed = false;
    }
    const bool spg_tail_allowed = !spg_post_case || handoff_cpp_tail_allow_spg_post_enabled();
    const bool allow_cpp_tail_reduction =
        cpp_tail_reduction_enabled() &&
        active_impl_mode() != GencanImplMode::kFortran &&
        !ab_compare_enabled() &&
        !tn_post_handoff_unsafe_enabled() &&
        !spg_post_handoff_enabled();
    const bool handoff_cpp_tail =
        handoff_cpp_tail_enabled() &&
        active_impl_mode() == GencanImplMode::kCpp &&
        !ab_compare_enabled() &&
        (!init1_active || handoff_cpp_tail_allow_init1_enabled()) &&
        tn_tail_allowed &&
        spg_tail_allowed &&
        (tn_handoff_requested || spg_post_handoff);
    const bool handoff_cpp_tail_init1_guarded =
        handoff_cpp_tail_enabled() &&
        active_impl_mode() == GencanImplMode::kCpp &&
        !ab_compare_enabled() &&
        init1_active &&
        !handoff_cpp_tail_allow_init1_enabled() &&
        (tn_handoff_requested || spg_post_handoff);
    const bool handoff_cpp_tail_spg_guarded =
        handoff_cpp_tail_enabled() &&
        active_impl_mode() == GencanImplMode::kCpp &&
        !ab_compare_enabled() &&
        spg_post_case &&
        !handoff_cpp_tail_allow_spg_post_enabled() &&
        (tn_handoff_requested || spg_post_handoff);
    const bool suppress_fallback_marker =
        (tn_post_case && tn_handoff_requested) ||
        spg_post_case ||
        allow_cpp_tail_reduction;
    const bool force_cpp_tn_post_tail =
        tn_post_case &&
        active_impl_mode() == GencanImplMode::kCpp &&
        !ab_compare_enabled() &&
        cpp_tail_reduction_enabled() &&
        !tn_handoff_requested &&
        !spg_post_handoff_enabled();
    const bool force_cpp_tail_in_cpp_mode =
        active_impl_mode() == GencanImplMode::kCpp &&
        !ab_compare_enabled() &&
        !tn_handoff_requested &&
        !spg_post_handoff_enabled();
    const bool force_cpp_tail_in_ab_init1_tn_post =
        active_impl_mode() != GencanImplMode::kFortran &&
        ab_compare_enabled() &&
        init1_active &&
        tn_post_case;
    const bool force_cpp_tail_in_ab_mode =
        active_impl_mode() == GencanImplMode::kAb &&
        ab_compare_enabled() &&
        !init1_active &&
        !tn_handoff_requested &&
        !spg_post_handoff_enabled();
    const bool effective_block_tail =
        block_cpp_tail_enabled() ||
        allow_cpp_tail_reduction ||
        handoff_cpp_tail ||
        force_cpp_tn_post_tail ||
        force_cpp_tail_in_ab_init1_tn_post ||
        force_cpp_tail_in_cpp_mode ||
        force_cpp_tail_in_ab_mode;
    return TailFallbackPolicy{
        tn_post_case,
        spg_post_case,
        tn_handoff_requested,
        spg_post_handoff,
        allow_cpp_tail_reduction,
        handoff_cpp_tail,
        handoff_cpp_tail_init1_guarded,
        handoff_cpp_tail_spg_guarded,
        suppress_fallback_marker,
        effective_block_tail
    };
}

struct FallbackPostDebugInfo {
    bool captured = false;
    int line_inform = 0;
    int post_inform = 0;
    int nind = 0;
    double f_before = 0.0;
    double f_after = 0.0;
    double gpsupn = 0.0;
    double gpeucn2 = 0.0;
};

void emit_fallback_debug_cpp(
    const TailFallbackPolicy& tail_policy,
    const char* fallback_reason,
    const FallbackPostDebugInfo& tn_post_debug,
    const FallbackPostDebugInfo& spg_post_debug
) {
    if (!gencan_debug_enabled()) {
        return;
    }
    std::fprintf(
        stderr,
        "[gencan-cpp-tail-policy] reason=%s mode=%d tn_case=%d spg_case=%d tn_handoff_req=%d spg_handoff=%d allow_cpp_tail_reduction=%d handoff_cpp_tail=%d effective_block_tail=%d\n",
        fallback_reason,
        static_cast<int>(active_impl_mode()),
        tail_policy.tn_post_case ? 1 : 0,
        tail_policy.spg_post_case ? 1 : 0,
        tail_policy.tn_handoff_requested ? 1 : 0,
        tail_policy.spg_post_handoff ? 1 : 0,
        tail_policy.allow_cpp_tail_reduction ? 1 : 0,
        tail_policy.handoff_cpp_tail ? 1 : 0,
        tail_policy.effective_block_tail ? 1 : 0
    );
    const bool tn_handoff_active = tail_policy.tn_post_case && tail_policy.tn_handoff_requested;
    const bool tn_unsafe = tn_post_handoff_unsafe_enabled();
    if (tn_handoff_active) {
        std::fprintf(stderr, "[gencan-cpp-handoff] reason=%s mode=%d\n", fallback_reason, static_cast<int>(active_impl_mode()));
        bool replay = tn_post_handoff_cpp_replay_enabled();
        if (tn_unsafe && !replay) {
            std::fprintf(stderr, "[gencan-cpp-handoff-canonicalize-cpp-replay-forced] reason=%s mode=%d\n", fallback_reason, static_cast<int>(active_impl_mode()));
            replay = true;
        }
        if (tn_post_handoff_canonicalize_enabled() && replay) {
            std::fprintf(stderr, "[gencan-cpp-handoff-canonicalize-cpp-replay] reason=%s mode=%d\n", fallback_reason, static_cast<int>(active_impl_mode()));
        }
    }
    if (tail_policy.spg_post_handoff) {
        std::fprintf(stderr, "[gencan-cpp-handoff] reason=spg_post_nonterminal mode=%d\n", static_cast<int>(active_impl_mode()));
        std::fprintf(stderr, "[gencan-cpp-handoff-cpp] reason=spg_post_nonterminal mode=%d\n", static_cast<int>(active_impl_mode()));
    }
    if (tail_policy.handoff_cpp_tail) {
        std::fprintf(stderr, "[gencan-cpp-handoff-cpp-tail] reason=%s mode=%d\n", fallback_reason, static_cast<int>(active_impl_mode()));
    }
    if (tail_policy.handoff_cpp_tail_init1_guarded) {
        std::fprintf(
            stderr,
            "[gencan-cpp-handoff-cpp-tail-init1-guard] reason=%s mode=%d\n",
            fallback_reason,
            static_cast<int>(active_impl_mode())
        );
    }
    if (tail_policy.handoff_cpp_tail_spg_guarded) {
        std::fprintf(
            stderr,
            "[gencan-cpp-handoff-cpp-tail-spg-guard] reason=%s mode=%d\n",
            fallback_reason,
            static_cast<int>(active_impl_mode())
        );
    }
    if (!tail_policy.suppress_fallback_marker) {
        std::fprintf(
            stderr,
            "[gencan-cpp-fallback] reason=%s mode=%d\n",
            fallback_reason,
            static_cast<int>(active_impl_mode())
        );
    }
    if (std::string(fallback_reason) == "tn_post_nonterminal" && tn_post_debug.captured) {
        std::fprintf(
            stderr,
            "[gencan-cpp-fallback-tn-post] line_inform=%d post_inform=%d nind=%d f_before=%.16e f_after=%.16e gpsupn=%.16e gpeucn2=%.16e\n",
            tn_post_debug.line_inform,
            tn_post_debug.post_inform,
            tn_post_debug.nind,
            tn_post_debug.f_before,
            tn_post_debug.f_after,
            tn_post_debug.gpsupn,
            tn_post_debug.gpeucn2
        );
    }
    if (std::string(fallback_reason) == "spg_post_nonterminal" && spg_post_debug.captured) {
        std::fprintf(
            stderr,
            "[gencan-cpp-fallback-spg-post] line_inform=%d post_inform=%d nind=%d f_before=%.16e f_after=%.16e gpsupn=%.16e gpeucn2=%.16e\n",
            spg_post_debug.line_inform,
            spg_post_debug.post_inform,
            spg_post_debug.nind,
            spg_post_debug.f_before,
            spg_post_debug.f_after,
            spg_post_debug.gpsupn,
            spg_post_debug.gpeucn2
        );
    }
}

void maybe_apply_fallback_seed_state_cpp(
    const bool fallback_seed_state,
    const bool fallback_x_seed_valid,
    const int n_val,
    const std::vector<double>& fallback_x_seed,
    const char* fallback_x_seed_reason,
    double* x
) {
    if (!fallback_seed_state || !fallback_x_seed_valid) {
        return;
    }
    for (int i = 0; i < n_val; ++i) {
        x[i] = fallback_x_seed[i];
    }
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp-fallback-seed] reason=%s mode=%d\n",
            fallback_x_seed_reason,
            static_cast<int>(active_impl_mode())
        );
    }
}

void maybe_apply_fallback_gradient_seed_state_cpp(
    const bool fallback_seed_state,
    const bool fallback_g_seed_valid,
    const int n_val,
    const std::vector<double>& fallback_g_seed,
    const char* fallback_g_seed_reason,
    double* g
) {
    if (!fallback_seed_state || !fallback_g_seed_valid) {
        return;
    }
    for (int i = 0; i < n_val; ++i) {
        g[i] = fallback_g_seed[i];
    }
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp-fallback-seed-g] reason=%s mode=%d\n",
            fallback_g_seed_reason,
            static_cast<int>(active_impl_mode())
        );
    }
}

void maybe_apply_fallback_f_seed_state_cpp(
    const bool fallback_f_seed_state,
    const bool fallback_f_seed_valid,
    const double fallback_f_seed,
    const char* fallback_f_seed_reason,
    double* f
) {
    if (!fallback_f_seed_state || !fallback_f_seed_valid) {
        return;
    }
    *f = fallback_f_seed;
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp-fallback-seed-f] reason=%s mode=%d f=%.16e\n",
            fallback_f_seed_reason,
            static_cast<int>(active_impl_mode()),
            fallback_f_seed
        );
    }
}

bool packmolprecision_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    int* eval_flag
) {
    double f_unused = 0.0;
    packmol_evalal_fortran_c(n, x, m, lambda, rho, &f_unused, eval_flag);
    if (*eval_flag < 0) {
        return false;
    }
    return packmolprecision_from_state_cpp();
}

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
) {
    const int nind_val = *nind;
    const int n_val = *n;

    *inform = 0;
    double xsupn = 0.0;
    double dsupn = 0.0;
    for (int i = 0; i < nind_val; ++i) {
        xsupn = std::max(xsupn, std::abs(x[i]));
        dsupn = std::max(dsupn, std::abs(d[i]));
    }
    if (dsupn < 1.0e-20) {
        dsupn = 1.0e-20;
    }
    const double step = std::max((*sterel) * xsupn, *steabs) / dsupn;

    for (int i = 0; i < nind_val; ++i) {
        xtmp[i] = x[i] + step * d[i];
    }

    calcg_cpp_reduced(nind, ind, xtmp, n, x, m, lambda, rho, gtype, hd, sterel, steabs, inform);
    if (*inform < 0) {
        return;
    }

    for (int i = 0; i < nind_val; ++i) {
        hd[i] = (hd[i] - g[i]) / step;
    }
    for (int i = nind_val; i < n_val; ++i) {
        hd[i] = 0.0;
    }
}

void calchd_cpp_reduced(
    const int* nind,
    const int* ind,
    double* x,
    double* d,
    double* g,
    const int* n,
    const double* xc,
    double* hd
) {
    const int nind_val = *nind;
    const int n_val = *n;

    for (int i = nind_val; i < n_val; ++i) {
        d[i] = 0.0;
        x[i] = xc[i];
    }

    expand_inplace(nind_val, ind, x);
    expand_inplace(nind_val, ind, d);
    expand_inplace(nind_val, ind, g);

    // Keep legacy calchd semantics: this wrapper performs reduced/full-space
    // packing and leaves hd values to external evalhd side effects.
    shrink_inplace(nind_val, ind, x);
    shrink_inplace(nind_val, ind, d);
    shrink_inplace(nind_val, ind, g);
    shrink_inplace(nind_val, ind, hd);
}

void spgls_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    const double* g,
    const double* l,
    const double* u,
    const double* lamspg,
    const double* nint,
    const int* mininterp,
    const double* fmin,
    const int* maxfc,
    int* fcnt,
    int* inform,
    double* xtrial,
    double* d,
    const double* gamma,
    const double* sigma1,
    const double* sigma2,
    const double* epsrel,
    const double* epsabs
) {
    const int n_val = *n;
    const int m_val = *m;
    const double f_ref = *f;
    const double lamspg_val = *lamspg;
    const double nint_val = *nint;
    const double sigma1_val = *sigma1;
    const double sigma2_val = *sigma2;
    const double gamma_val = *gamma;
    const double epsrel_val = *epsrel;
    const double epsabs_val = *epsabs;
    const double fmin_val = *fmin;

    int interp = 0;
    double alpha = 1.0;
    for (int i = 0; i < n_val; ++i) {
        xtrial[i] = std::min(u[i], std::max(l[i], x[i] - lamspg_val * g[i]));
        d[i] = xtrial[i] - x[i];
    }
    double gtd = 0.0;
    for (int i = 0; i < n_val; ++i) {
        gtd += g[i] * d[i];
    }

    double ftrial = 0.0;
    packmol_evalal_fortran_c(n, xtrial, m, lambda, rho, &ftrial, inform);
    *fcnt += 1;
    if (*inform < 0) {
        return;
    }

    while (true) {
        if (ftrial <= f_ref + gamma_val * alpha * gtd) {
            *f = ftrial;
            vec_copy(n_val, xtrial, x);
            *inform = 0;
            return;
        }

        if (ftrial <= fmin_val) {
            *f = ftrial;
            vec_copy(n_val, xtrial, x);
            *inform = 4;
            return;
        }

        if (*fcnt >= *maxfc) {
            if (ftrial < f_ref) {
                *f = ftrial;
                vec_copy(n_val, xtrial, x);
            }
            *inform = 8;
            return;
        }

        interp += 1;
        if (alpha < sigma1_val) {
            alpha = alpha / nint_val;
        } else {
            const double den = 2.0 * (ftrial - f_ref - alpha * gtd);
            double atmp = alpha / nint_val;
            if (den != 0.0) {
                atmp = (-gtd * alpha * alpha) / den;
            }
            if (atmp < sigma1_val || atmp > sigma2_val * alpha) {
                alpha = alpha / nint_val;
            } else {
                alpha = atmp;
            }
        }

        vec_trial_point(n_val, x, alpha, d, xtrial);

        packmol_evalal_fortran_c(&n_val, xtrial, &m_val, lambda, rho, &ftrial, inform);
        *fcnt += 1;
        if (*inform < 0) {
            return;
        }

        const bool samep = same_step_relative(n_val, alpha, d, x, epsrel_val, epsabs_val);

        if (interp >= *mininterp && samep) {
            if (ftrial < f_ref) {
                *f = ftrial;
                vec_copy(n_val, xtrial, x);
            }
            *inform = 6;
            return;
        }
    }
}

bool tnls_cpp_subset(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* l,
    const double* u,
    double* f,
    double* g,
    const double* d,
    const double* amax,
    const int* rbdtype,
    const int* rbdind,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    const int* maxfc,
    const int* gtype,
    int* fcnt,
    int* gcnt,
    int* intcnt,
    int* exgcnt,
    int* exbcnt,
    int* inform,
    double* xplus,
    double* xtmp,
    double* xbext,
    const double* gamma,
    const double* beta,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs
) {
    const bool valid_bound =
        (*rbdtype == 1 || *rbdtype == 2) && (*rbdind >= 1 && *rbdind <= *nind);
    const bool has_active_bound = valid_bound;

    const int nind_val = *nind;
    const int n_val = *n;
    const int m_val = *m;
    const double f0 = *f;

    auto compute_gradient = [&]() -> bool {
        calcg_cpp_reduced(
            nind, ind, x, &n_val, x, &m_val, lambda, rho, gtype, g, sterel, steabs, inform
        );
        *gcnt += 1;
        return *inform >= 0;
    };

    auto finish_extrapolation = [&](double fbext, int code, bool need_gradient_update) -> bool {
        if (need_gradient_update) {
            if (compute_gradient()) {
                if (*f < fbext) {
                    *exgcnt += 1;
                } else {
                    *exbcnt += 1;
                }
            }
        }
        *inform = code;
        return true;
    };

    auto run_extrapolation = [&](double& alpha_ref, double& fplus_ref, bool gradient_already_at_xplus) -> bool {
        double fbext = fplus_ref;
        for (int i = 0; i < nind_val; ++i) {
            xbext[i] = xplus[i];
        }
        int extrap = 0;

        while (true) {
            const bool need_gradient_update = (extrap != 0) || (*amax <= 1.0) || (!gradient_already_at_xplus);

            if (fplus_ref <= *fmin) {
                *f = fplus_ref;
                vec_copy(nind_val, xplus, x);
                return finish_extrapolation(fbext, 4, need_gradient_update);
            }

            if (*fcnt >= *maxfc) {
                *f = fplus_ref;
                vec_copy(nind_val, xplus, x);
                return finish_extrapolation(fbext, 8, need_gradient_update);
            }

            if (extrap >= *maxextrap) {
                *f = fplus_ref;
                vec_copy(nind_val, xplus, x);
                return finish_extrapolation(fbext, 7, need_gradient_update);
            }

            double atmp = (*next) * alpha_ref;
            if (alpha_ref < *amax && atmp > *amax) {
                atmp = *amax;
            }

            vec_trial_point(nind_val, x, atmp, d, xtmp);
            if (atmp == *amax && has_active_bound) {
                const int idx = *rbdind - 1;
                xtmp[idx] = (*rbdtype == 1) ? l[idx] : u[idx];
            }
            if (atmp > *amax) {
                for (int i = 0; i < nind_val; ++i) {
                    xtmp[i] = std::max(l[i], std::min(xtmp[i], u[i]));
                }
            }

            if (alpha_ref > *amax) {
                const bool samep =
                    same_point_relative(nind_val, xtmp, xplus, xplus, *epsrel, *epsabs);
                if (samep) {
                    *f = fplus_ref;
                    vec_copy(nind_val, xplus, x);
                    return finish_extrapolation(fbext, 0, need_gradient_update);
                }
            }

            double ftmp = 0.0;
            calcf_cpp_reduced(nind, ind, xtmp, &n_val, x, &m_val, lambda, rho, &ftmp, inform);
            *fcnt += 1;
            if (*inform < 0) {
                *f = fbext;
                vec_copy(nind_val, xbext, x);
                if (need_gradient_update && compute_gradient()) {
                    *exbcnt += 1;
                }
                *inform = 0;
                return true;
            }

            if (ftmp < fplus_ref) {
                alpha_ref = atmp;
                fplus_ref = ftmp;
                vec_copy(nind_val, xtmp, xplus);
                extrap += 1;
            } else {
                *f = fplus_ref;
                vec_copy(nind_val, xplus, x);
                return finish_extrapolation(fbext, 0, need_gradient_update);
            }
        }
        };

    double gtd = 0.0;
    for (int i = 0; i < nind_val; ++i) {
        gtd += g[i] * d[i];
    }

    double alpha = std::min(1.0, *amax);
    vec_trial_point(nind_val, x, alpha, d, xplus);
    if (alpha == *amax && has_active_bound) {
        const int idx = *rbdind - 1;
        xplus[idx] = (*rbdtype == 1) ? l[idx] : u[idx];
    }

    double fplus = 0.0;
    calcf_cpp_reduced(nind, ind, xplus, n, x, m, lambda, rho, &fplus, inform);
    *fcnt += 1;
    if (*inform < 0) {
        return true;
    }

    if (*amax > 1.0) {
        if (fplus <= *f + *gamma * alpha * gtd) {
            calcg_cpp_reduced(
                nind, ind, xplus, &n_val, x, &m_val, lambda, rho, gtype, g, sterel, steabs, inform
            );
            *gcnt += 1;
            if (*inform < 0) {
                return true;
            }

            double gptd = 0.0;
            for (int i = 0; i < nind_val; ++i) {
                gptd += g[i] * d[i];
            }

            if (gptd < (*beta) * gtd) {
                return run_extrapolation(alpha, fplus, true);
            }

            *f = fplus;
            vec_copy(nind_val, xplus, x);
            *inform = 0;
            return true;
        }
        // Armijo does not hold at the interior trial.
        // Continue with interpolation logic below in C++.
    }

    if (fplus < *f) {
        return run_extrapolation(alpha, fplus, false);
    }

    *intcnt += 1;
    int interp = 0;
    while (true) {
        if (fplus <= *fmin) {
            *f = fplus;
            vec_copy(nind_val, xplus, x);
            if (!compute_gradient()) {
                return true;
            }
            *inform = 4;
            return true;
        }

        if (*fcnt >= *maxfc) {
            if (fplus < f0) {
                *f = fplus;
                vec_copy(nind_val, xplus, x);
                if (!compute_gradient()) {
                    return true;
                }
            }
            *inform = 8;
            return true;
        }

        if (fplus <= *f + *gamma * alpha * gtd) {
            *f = fplus;
            vec_copy(nind_val, xplus, x);
            if (!compute_gradient()) {
                return true;
            }
            *inform = 0;
            return true;
        }

        interp += 1;
        if (alpha < *sigma1) {
            alpha = alpha / *nint;
        } else {
            const double den = 2.0 * (fplus - *f - alpha * gtd);
            double atmp = alpha / *nint;
            if (den != 0.0) {
                atmp = (-gtd * alpha * alpha) / den;
            }
            if (atmp < *sigma1 || atmp > *sigma2 * alpha) {
                alpha = alpha / *nint;
            } else {
                alpha = atmp;
            }
        }

        vec_trial_point(nind_val, x, alpha, d, xplus);
        calcf_cpp_reduced(nind, ind, xplus, &n_val, x, &m_val, lambda, rho, &fplus, inform);
        *fcnt += 1;
        if (*inform < 0) {
            return true;
        }

        const bool samep = same_step_relative(nind_val, alpha, d, x, *epsrel, *epsabs);
        if (interp >= *mininterp && samep) {
            *inform = 6;
            return true;
        }
    }
}

void cg_cpp_full(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* eps,
    const double* epsnqmp,
    const int* maxitnqmp,
    const int* maxit,
    const bool* nearlyq,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const double* theta,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infabs,
    double* s,
    double* w,
    double* y,
    double* r,
    double* d,
    double* sprev,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform
) {
    const int nind_val = *nind;
    (void)iprint;
    const bool debug = gencan_debug_enabled();
    *rbdtype = 0;
    *rbdind = 0;
    auto norm2 = [&](const double* v) { return norm2_kernel(nind_val, v); };
    auto dot = [&](const double* a, const double* b) { return dot_kernel(nind_val, a, b); };

    const double gnorm2 = norm2(g);

    *iter = 0;
    int itnqmp = 0;
    double qprev = *infabs;
    double bestprog = 0.0;

    for (int i = 0; i < nind_val; ++i) {
        s[i] = 0.0;
        r[i] = g[i];
    }

    double q = 0.0;
    double snorm2 = 0.0;
    double rnorm2 = gnorm2;
    double dnorm2 = 0.0;
    double dtr = 0.0;
    double alpha = 0.0;
    double dtw = 0.0;
    double rnorm2prev = rnorm2;

    while (true) {
        if (rnorm2 <= 1.0e-16 ||
            (((rnorm2 <= (*eps) * (*eps) * gnorm2) || (rnorm2 <= 1.0e-10 && *iter != 0)) &&
             *iter >= 4)) {
            *rbdtype = 0;
            *rbdind = 0;
            *inform = 0;
            return;
        }

        if (*iter >= std::max(4, *maxit)) {
            *rbdtype = 0;
            *rbdind = 0;
            *inform = 8;
            return;
        }

        if (*iter == 0) {
            for (int i = 0; i < nind_val; ++i) {
                d[i] = -r[i];
            }
            dnorm2 = rnorm2;
            dtr = -rnorm2;
        } else {
            const double beta = rnorm2 / rnorm2prev;
            for (int i = 0; i < nind_val; ++i) {
                d[i] = -r[i] + beta * d[i];
            }
            dnorm2 = rnorm2 - 2.0 * beta * (dtr + alpha * dtw) + beta * beta * dnorm2;
            dtr = -rnorm2 + beta * (dtr + alpha * dtw);
        }

        if (dtr > 0.0) {
            for (int i = 0; i < nind_val; ++i) {
                d[i] = -d[i];
            }
            dtr = -dtr;
        }

        if (*htvtype == 0) {
            calchd_cpp_reduced(nind, ind, x, d, g, n, x, w);
        } else if (*htvtype == 1) {
            calchddiff_cpp_reduced(
                nind, ind, n, x, d, g, m, lambda, rho, gtype, w, y, sterel, steabs, inform
            );
        } else {
            *inform = -1;
            return;
        }

        if (*inform < 0) {
            return;
        }

        dtw = dot(d, w);

        if (debug) {
            std::fprintf(
                stderr,
                "[gencan-cg-iter] call=%d iter=%d rnorm2=%.16e dnorm2=%.16e dtr=%.16e dtw=%.16e\n",
                g_cg_active_call_id, *iter, rnorm2, dnorm2, dtr, dtw
            );
        }

        const double dts = dot(d, s);

        double amax1 = *infabs;
        double amax1n = -*infabs;
        if (*trtype == 0) {
            const double aa = dnorm2;
            const double bb = 2.0 * dts;
            const double cc = snorm2 - (*delta) * (*delta);
            const double dd = std::sqrt(bb * bb - 4.0 * aa * cc);
            amax1 = (-bb + dd) / (2.0 * aa);
            amax1n = (-bb - dd) / (2.0 * aa);
        } else if (*trtype == 1) {
            for (int i = 0; i < nind_val; ++i) {
                if (d[i] > 0.0) {
                    amax1 = std::min(amax1, ((*delta) - s[i]) / d[i]);
                    amax1n = std::max(amax1n, ((-*delta) - s[i]) / d[i]);
                } else if (d[i] < 0.0) {
                    amax1 = std::min(amax1, ((-*delta) - s[i]) / d[i]);
                    amax1n = std::max(amax1n, ((*delta) - s[i]) / d[i]);
                }
            }
        }

        double amax2 = *infabs;
        double amax2n = -*infabs;
        int rbdposaind = 0;
        int rbdposatype = 0;
        int rbdnegaind = 0;
        int rbdnegatype = 0;
        for (int i = 0; i < nind_val; ++i) {
            if (d[i] > 0.0) {
                const double amax2x = (u[i] - x[i] - s[i]) / d[i];
                if (amax2x < amax2) {
                    amax2 = amax2x;
                    rbdposaind = i + 1;
                    rbdposatype = 2;
                }
                const double amax2nx = (l[i] - x[i] - s[i]) / d[i];
                if (amax2nx > amax2n) {
                    amax2n = amax2nx;
                    rbdnegaind = i + 1;
                    rbdnegatype = 1;
                }
            } else if (d[i] < 0.0) {
                const double amax2x = (l[i] - x[i] - s[i]) / d[i];
                if (amax2x < amax2) {
                    amax2 = amax2x;
                    rbdposaind = i + 1;
                    rbdposatype = 1;
                }
                const double amax2nx = (u[i] - x[i] - s[i]) / d[i];
                if (amax2nx > amax2n) {
                    amax2n = amax2nx;
                    rbdnegaind = i + 1;
                    rbdnegatype = 2;
                }
            }
        }

        const double amax = std::min(amax1, amax2);
        const double amaxn = std::max(amax1n, amax2n);

        qprev = q;
        if (dtw > 0.0) {
            alpha = std::min(amax, rnorm2 / dtw);
            q = q + 0.5 * alpha * alpha * dtw + alpha * dtr;
        } else {
            const double qamax = q + 0.5 * amax * amax * dtw + amax * dtr;
            if (*iter == 0) {
                alpha = amax;
                q = qamax;
            } else {
                const double qamaxn = q + 0.5 * amaxn * amaxn * dtw + amaxn * dtr;
                const bool allow_negative_curvature_step = *nearlyq || cg_dtw_relax_enabled();
                if (allow_negative_curvature_step && (qamax < q || qamaxn < q)) {
                    if (qamax < qamaxn) {
                        alpha = amax;
                        q = qamax;
                    } else {
                        alpha = amaxn;
                        q = qamaxn;
                    }
                } else {
                    *rbdtype = 0;
                    *rbdind = 0;
                    *inform = 7;
                    if (debug) {
                        std::fprintf(
                            stderr,
                            "[gencan-cg-exit] call=%d iter=%d reason=inform7 dtw=%.16e amax=%.16e amaxn=%.16e q=%.16e qamax=%.16e qamaxn=%.16e\n",
                            g_cg_active_call_id, *iter, dtw, amax, amaxn, q, qamax, qamaxn
                        );
                    }
                    return;
                }
            }
        }

        if (debug) {
            std::fprintf(
                stderr,
                "[gencan-cg-step] call=%d iter=%d alpha=%.16e amax=%.16e amax1=%.16e amax2=%.16e amaxn=%.16e q=%.16e\n",
                g_cg_active_call_id, *iter, alpha, amax, amax1, amax2, amaxn, q
            );
        }

        for (int i = 0; i < nind_val; ++i) {
            sprev[i] = s[i];
            s[i] = s[i] + alpha * d[i];
        }
        const double snorm2prev = snorm2;
        snorm2 = snorm2 + alpha * alpha * dnorm2 + 2.0 * alpha * dts;

        rnorm2prev = rnorm2;
        for (int i = 0; i < nind_val; ++i) {
            r[i] = r[i] + alpha * w[i];
        }
        rnorm2 = norm2(r);

        *iter += 1;

        const double gts = dot(g, s);
        if (gts > 0.0 || gts * gts < (*theta) * (*theta) * gnorm2 * snorm2) {
            for (int i = 0; i < nind_val; ++i) {
                s[i] = sprev[i];
            }
            snorm2 = snorm2prev;
            q = qprev;
            *rbdtype = 0;
            *rbdind = 0;
            *inform = 3;
            if (debug) {
                std::fprintf(
                    stderr,
                    "[gencan-cg-exit] call=%d iter=%d reason=angle gts=%.16e gnorm2=%.16e snorm2=%.16e\n",
                    g_cg_active_call_id, *iter, gts, gnorm2, snorm2
                );
            }
            return;
        }

        if (alpha == amax2 || alpha == amax2n) {
            if (alpha == amax2) {
                *rbdind = rbdposaind;
                *rbdtype = rbdposatype;
            } else {
                *rbdind = rbdnegaind;
                *rbdtype = rbdnegatype;
            }
            *inform = 2;
            if (debug) {
                std::fprintf(
                    stderr,
                    "[gencan-cg-exit] call=%d iter=%d reason=box alpha=%.16e amax2=%.16e amax2n=%.16e rbdtype=%d rbdind=%d\n",
                    g_cg_active_call_id, *iter, alpha, amax2, amax2n, *rbdtype, *rbdind
                );
            }
            return;
        }

        if (alpha == amax1 || alpha == amax1n) {
            *rbdtype = 0;
            *rbdind = 0;
            *inform = 1;
            if (debug) {
                std::fprintf(
                    stderr,
                    "[gencan-cg-exit] call=%d iter=%d reason=trust alpha=%.16e amax1=%.16e amax1n=%.16e\n",
                    g_cg_active_call_id, *iter, alpha, amax1, amax1n
                );
            }
            return;
        }

        bool samep = true;
        for (int i = 0; i < nind_val; ++i) {
            if (std::abs(alpha * d[i]) > std::max((*epsrel) * std::abs(s[i]), *epsabs)) {
                samep = false;
            }
        }
        if (samep) {
            *rbdtype = 0;
            *rbdind = 0;
            *inform = 6;
            if (debug) {
                std::fprintf(
                    stderr,
                    "[gencan-cg-exit] call=%d iter=%d reason=samep alpha=%.16e\n",
                    g_cg_active_call_id, *iter, alpha
                );
            }
            return;
        }

        const double currprog = qprev - q;
        bestprog = std::max(currprog, bestprog);
        if (currprog <= (*epsnqmp) * bestprog) {
            itnqmp += 1;
            if (itnqmp >= *maxitnqmp) {
                *rbdtype = 0;
                *rbdind = 0;
                *inform = 4;
                if (debug) {
                    std::fprintf(
                        stderr,
                        "[gencan-cg-exit] call=%d iter=%d reason=nqmp currprog=%.16e bestprog=%.16e epsnqmp=%.16e itnqmp=%d\n",
                        g_cg_active_call_id, *iter, currprog, bestprog, *epsnqmp, itnqmp
                    );
                }
                return;
            }
        } else {
            itnqmp = 0;
        }
    }
}

}  // namespace

extern "C" void packmol_gencan_spgls_bridge(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    const double* g,
    const double* l,
    const double* u,
    const double* lamspg,
    const double* nint,
    const int* mininterp,
    const double* fmin,
    const int* maxfc,
    const int* iprint,
    int* fcnt,
    int* inform,
    double* xtrial,
    double* d,
    const double* gamma,
    const double* sigma1,
    const double* sigma2,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
) {
    // Phase I runtime switch scaffold:
    // - fortran: force legacy kernel
    // - cpp: replacement path (currently forwards to legacy kernel)
    // - ab: A/B harness mode (candidate currently legacy kernel placeholder)
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_spgls_fortran_c(
                n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp,
                fmin, maxfc, iprint, fcnt, inform, xtrial, d, gamma, sigma1, sigma2,
                sterel, steabs, epsrel, epsabs, infrel, infabs
            );
            return;
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            // First real replacement kernel in C++.
            // This preserves legacy callback behavior through packmol_evalal_fortran_c.
            spgls_cpp(
                n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, fmin,
                maxfc, fcnt, inform, xtrial, d, gamma, sigma1, sigma2, epsrel, epsabs
            );
            return;
    }
}

extern "C" void packmol_tnls_fortran_c(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* l,
    const double* u,
    double* f,
    double* g,
    const double* d,
    const double* amax,
    const int* rbdtype,
    const int* rbdind,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    const int* maxfc,
    const int* gtype,
    const int* iprint,
    int* fcnt,
    int* gcnt,
    int* intcnt,
    int* exgcnt,
    int* exbcnt,
    int* inform,
    double* xplus,
    double* xtmp,
    double* xbext,
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

extern "C" void packmol_gencan_tnls_bridge(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    const double* l,
    const double* u,
    double* f,
    double* g,
    const double* d,
    const double* amax,
    const int* rbdtype,
    const int* rbdind,
    const double* nint,
    const double* next,
    const int* mininterp,
    const int* maxextrap,
    const double* fmin,
    const int* maxfc,
    const int* gtype,
    const int* iprint,
    int* fcnt,
    int* gcnt,
    int* intcnt,
    int* exgcnt,
    int* exbcnt,
    int* inform,
    double* xplus,
    double* xtmp,
    double* xbext,
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
) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_tnls_fortran_c(
                nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, rbdind,
                nint, next, mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, gcnt,
                intcnt, exgcnt, exbcnt, inform, xplus, xtmp, xbext, gamma, beta, sigma1,
                sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
            );
            return;
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            if (tnls_cpp_subset(
                    nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, rbdind,
                    nint, next, mininterp, maxextrap, fmin, maxfc, gtype, fcnt, gcnt, intcnt,
                    exgcnt, exbcnt, inform, xplus, xtmp, xbext, gamma, beta, sigma1, sigma2,
                    sterel, steabs, epsrel, epsabs)) {
                return;
            }
            return;
    }
}

extern "C" void packmol_cg_fortran_c(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* eps,
    const double* epsnqmp,
    const int* maxitnqmp,
    const int* maxit,
    const bool* nearlyq,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* s,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform,
    double* w,
    double* y,
    double* r,
    double* d,
    double* sprev,
    const double* theta,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
);

extern "C" void packmol_gencan_cg_bridge(
    const int* nind,
    const int* ind,
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    const double* delta,
    const double* l,
    const double* u,
    const double* eps,
    const double* epsnqmp,
    const int* maxitnqmp,
    const int* maxit,
    const bool* nearlyq,
    const int* gtype,
    const int* htvtype,
    const int* trtype,
    const int* iprint,
    const int* ncomp,
    double* s,
    int* iter,
    int* rbdtype,
    int* rbdind,
    int* inform,
    double* w,
    double* y,
    double* r,
    double* d,
    double* sprev,
    const double* theta,
    const double* sterel,
    const double* steabs,
    const double* epsrel,
    const double* epsabs,
    const double* infrel,
    const double* infabs
) {
    static int cg_call_counter = 0;
    cg_call_counter += 1;
    if (gencan_debug_enabled()) {
        double gnorm2_in = 0.0;
        double xsum_in = 0.0;
        for (int i = 0; i < *nind; ++i) {
            gnorm2_in += g[i] * g[i];
            xsum_in += x[i];
        }
        std::fprintf(
            stderr,
            "[gencan-cg-entry] call=%d mode=%d nind=%d gnorm2=%.16e xsum=%.16e delta=%.16e eps=%.16e epsnqmp=%.16e maxit=%d maxitnqmp=%d\n",
            cg_call_counter, static_cast<int>(active_impl_mode()), *nind, gnorm2_in, xsum_in,
            *delta, *eps, *epsnqmp, *maxit, *maxitnqmp
        );
    }
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_cg_fortran_c(
                nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp,
                maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, ncomp, s,
                iter, rbdtype, rbdind, inform, w, y, r, d, sprev, theta, sterel,
                steabs, epsrel, epsabs, infrel, infabs
            );
            if (gencan_debug_enabled()) {
                double snorm2 = 0.0;
                for (int i = 0; i < *nind; ++i) {
                    snorm2 += s[i] * s[i];
                }
                std::fprintf(
                    stderr,
                    "[gencan-cg-fortran] nind=%d nearlyq=%d iter=%d inform=%d rbdtype=%d rbdind=%d snorm2=%.16e\n",
                    *nind, *nearlyq ? 1 : 0, *iter, *inform, *rbdtype, *rbdind, snorm2
                );
            }
            return;
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb:
            int fortran_iter = 0;
            int fortran_rbdtype = 0;
            int fortran_rbdind = 0;
            int fortran_inform = 0;
            double fortran_snorm2 = 0.0;
            std::vector<double> s_ref_dbg;
            std::vector<double> r_ref_dbg;
            std::vector<double> d_ref_dbg;
            const bool shadow_enabled = cg_shadow_compare_enabled();
            std::vector<int> ind_ref;
            std::vector<double> x_ref, g_ref, s_ref, w_ref, y_ref, r_ref, d_ref, sprev_ref;
            if (shadow_enabled) {
                ind_ref.resize(*nind);
                x_ref.resize(*n);
                g_ref.resize(*n);
                s_ref.resize(*n);
                w_ref.resize(*n);
                y_ref.resize(*n);
                r_ref.resize(*n);
                d_ref.resize(*n);
                sprev_ref.resize(*n);
                for (int i = 0; i < *nind; ++i) ind_ref[i] = ind[i];
                for (int i = 0; i < *n; ++i) {
                    x_ref[i] = x[i];
                    g_ref[i] = g[i];
                    s_ref[i] = s[i];
                    w_ref[i] = w[i];
                    y_ref[i] = y[i];
                    r_ref[i] = r[i];
                    d_ref[i] = d[i];
                    sprev_ref[i] = sprev[i];
                }
            }
            if (shadow_enabled) {
                int iprint_shadow = *iprint;
                if (gencan_debug_enabled()) {
                    iprint_shadow = std::max(iprint_shadow, 4);
                }
                packmol_cg_fortran_c(
                    nind, ind_ref.data(), n, x_ref.data(), m, lambda, rho, g_ref.data(), delta, l, u, eps, epsnqmp,
                    maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, &iprint_shadow, ncomp, s_ref.data(),
                    &fortran_iter, &fortran_rbdtype, &fortran_rbdind, &fortran_inform, w_ref.data(), y_ref.data(),
                    r_ref.data(), d_ref.data(), sprev_ref.data(), theta, sterel, steabs, epsrel, epsabs, infrel, infabs
                );
                for (int i = 0; i < *nind; ++i) {
                    fortran_snorm2 += s_ref[i] * s_ref[i];
                }
                s_ref_dbg.assign(s_ref.begin(), s_ref.end());
                r_ref_dbg.assign(r_ref.begin(), r_ref.end());
                d_ref_dbg.assign(d_ref.begin(), d_ref.end());
            }
            if (gencan_debug_enabled()) {
                std::fprintf(stderr, "[gencan-cg-cpp] nind=%d n=%d trtype=%d htvtype=%d gtype=%d\n",
                             *nind, *n, *trtype, *htvtype, *gtype);
            }
            g_cg_active_call_id = cg_call_counter;
            cg_cpp_full(
                nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp,
                maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, theta,
                sterel, steabs, epsrel, epsabs, infabs, s, w, y, r, d, sprev, iter,
                rbdtype, rbdind, inform
            );
            g_cg_active_call_id = 0;
            if (gencan_debug_enabled()) {
                double snorm2 = 0.0;
                for (int i = 0; i < *nind; ++i) {
                    snorm2 += s[i] * s[i];
                }
                std::fprintf(
                    stderr,
                    "[gencan-cg-cpp-out] nind=%d nearlyq=%d iter=%d inform=%d rbdtype=%d rbdind=%d snorm2=%.16e\n",
                    *nind, *nearlyq ? 1 : 0, *iter, *inform, *rbdtype, *rbdind, snorm2
                );
                if (shadow_enabled) {
                    double max_s_diff = 0.0;
                    double max_r_diff = 0.0;
                    double max_d_diff = 0.0;
                    int max_s_idx = -1;
                    int max_r_idx = -1;
                    int max_d_idx = -1;
                    for (int i = 0; i < *nind; ++i) {
                        if (i < static_cast<int>(s_ref_dbg.size())) {
                            const double ds = std::abs(s_ref_dbg[i] - s[i]);
                            if (ds > max_s_diff) {
                                max_s_diff = ds;
                                max_s_idx = i + 1;
                            }
                        }
                        if (i < static_cast<int>(r_ref_dbg.size())) {
                            const double dr = std::abs(r_ref_dbg[i] - r[i]);
                            if (dr > max_r_diff) {
                                max_r_diff = dr;
                                max_r_idx = i + 1;
                            }
                        }
                        if (i < static_cast<int>(d_ref_dbg.size())) {
                            const double dd = std::abs(d_ref_dbg[i] - d[i]);
                            if (dd > max_d_diff) {
                                max_d_diff = dd;
                                max_d_idx = i + 1;
                            }
                        }
                    }
                    std::fprintf(
                        stderr,
                        "[gencan-cg-shadow] nind=%d f_iter=%d c_iter=%d f_inform=%d c_inform=%d f_rbdtype=%d c_rbdtype=%d f_rbdind=%d c_rbdind=%d f_snorm2=%.16e c_snorm2=%.16e max|ds|=%.16e@%d max|dr|=%.16e@%d max|dd|=%.16e@%d\n",
                        *nind, fortran_iter, *iter, fortran_inform, *inform, fortran_rbdtype,
                        *rbdtype, fortran_rbdind, *rbdind, fortran_snorm2, snorm2,
                        max_s_diff, max_s_idx, max_r_diff, max_r_idx, max_d_diff, max_d_idx
                    );
                }
            }
            return;
    }
}

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

extern "C" void packmol_gencan_fortran_c(
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

extern "C" void packmol_gencan_gencan_bridge(
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

    packmol_gencan_gencan_bridge(
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

}  // namespace

extern "C" void packmol_gencan_easy_bridge(
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
        case GencanImplMode::kAb:
            easygencan_cpp(
                n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, iprint, ncomp,
                f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin
            );
            return;
    }
}

extern "C" void packmol_gencan_gencan_bridge(
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
) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_gencan_fortran_c(
                n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp,
                maxitngp, fmin, maxit, maxfc, udelta0, ucgmaxit, cgscre, cggpnf,
                cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq, nint, next, mininterp,
                maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2,
                gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt,
                tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe,
                inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi,
                theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs,
                infrel, infabs
            );
            return;
        case GencanImplMode::kCpp:
        case GencanImplMode::kAb: {
            const int n_val = *n;
            // Legacy Fortran treats these as output counters; caller values are
            // not part of the contract and may be uninitialized.
            *iter = 0;
            *fcnt = 0;
            *gcnt = 0;
            *cgcnt = 0;
            *spgiter = 0;
            *spgfcnt = 0;
            *tniter = 0;
            *tnfcnt = 0;
            *tnstpcnt = 0;
            *tnintcnt = 0;
            *tnexgcnt = 0;
            *tnexbcnt = 0;
            *tnintfe = 0;
            *tnexgfe = 0;
            *tnexbfe = 0;
            const char* fallback_reason = nullptr;
            std::vector<double> fallback_x_seed;
            bool fallback_x_seed_valid = false;
            const char* fallback_x_seed_reason = "";
            std::vector<double> fallback_g_seed;
            bool fallback_g_seed_valid = false;
            const char* fallback_g_seed_reason = "";
            double fallback_f_seed = 0.0;
            bool fallback_f_seed_valid = false;
            const char* fallback_f_seed_reason = "";
            bool fallback_cg_schedule_seed_valid = false;
            double fallback_cg_schedule_acgeps = 0.0;
            double fallback_cg_schedule_bcgeps = 0.0;
            double fallback_cg_schedule_gpeucn20 = 0.0;
            double fallback_cg_schedule_gpsupn0 = 0.0;
            bool gencan_cg_schedule_seed_valid = false;
            double gencan_cg_schedule_acgeps = 0.0;
            double gencan_cg_schedule_bcgeps = 0.0;
            double gencan_cg_schedule_gpeucn20 = 0.0;
            double gencan_cg_schedule_gpsupn0 = 0.0;
            bool fallback_counter_seed_valid = false;
            int fallback_iter_seed = 0;
            int fallback_fcnt_seed = 0;
            int fallback_gcnt_seed = 0;
            int fallback_cgcnt_seed = 0;
            int fallback_spgiter_seed = 0;
            int fallback_spgfcnt_seed = 0;
            int fallback_tniter_seed = 0;
            int fallback_tnfcnt_seed = 0;
            int fallback_tnstpcnt_seed = 0;
            int fallback_tnintcnt_seed = 0;
            int fallback_tnexgcnt_seed = 0;
            int fallback_tnexbcnt_seed = 0;
            int fallback_tnintfe_seed = 0;
            int fallback_tnexgfe_seed = 0;
            int fallback_tnexbfe_seed = 0;
            bool fallback_step_seed_valid = false;
            std::vector<double> fallback_s_seed;
            std::vector<double> fallback_y_seed;
            bool fallback_progress_seed_valid = false;
            double fallback_progress_fprev = 0.0;
            double fallback_progress_bestprog = 0.0;
            int fallback_progress_itnfp = 0;

            auto seed_fallback_state = [&](const char* reason, const std::vector<double>& x_seed,
                                           const std::vector<double>* g_seed, const double f_seed) {
                fallback_x_seed = x_seed;
                fallback_x_seed_valid = true;
                fallback_x_seed_reason = reason;
                if (g_seed != nullptr) {
                    fallback_g_seed = *g_seed;
                    fallback_g_seed_valid = true;
                    fallback_g_seed_reason = reason;
                } else {
                    fallback_g_seed.clear();
                    fallback_g_seed_valid = false;
                    fallback_g_seed_reason = "";
                }
                fallback_f_seed = f_seed;
                fallback_f_seed_valid = true;
                fallback_f_seed_reason = reason;
            };

            auto seed_fallback_counters = [&](const int iter_seed, const int fcnt_seed, const int gcnt_seed,
                                              const int cgcnt_seed, const int spgiter_seed,
                                              const int spgfcnt_seed, const int tniter_seed,
                                              const int tnfcnt_seed, const int tnstpcnt_seed,
                                              const int tnintcnt_seed, const int tnexgcnt_seed,
                                              const int tnexbcnt_seed, const int tnintfe_seed,
                                              const int tnexgfe_seed, const int tnexbfe_seed) {
                fallback_counter_seed_valid = true;
                fallback_iter_seed = iter_seed;
                fallback_fcnt_seed = fcnt_seed;
                fallback_gcnt_seed = gcnt_seed;
                fallback_cgcnt_seed = cgcnt_seed;
                fallback_spgiter_seed = spgiter_seed;
                fallback_spgfcnt_seed = spgfcnt_seed;
                fallback_tniter_seed = tniter_seed;
                fallback_tnfcnt_seed = tnfcnt_seed;
                fallback_tnstpcnt_seed = tnstpcnt_seed;
                fallback_tnintcnt_seed = tnintcnt_seed;
                fallback_tnexgcnt_seed = tnexgcnt_seed;
                fallback_tnexbcnt_seed = tnexbcnt_seed;
                fallback_tnintfe_seed = tnintfe_seed;
                fallback_tnexgfe_seed = tnexgfe_seed;
                fallback_tnexbfe_seed = tnexbfe_seed;
            };

            auto seed_fallback_steps = [&](const std::vector<double>& s_seed, const std::vector<double>& y_seed) {
                fallback_step_seed_valid = true;
                fallback_s_seed = s_seed;
                fallback_y_seed = y_seed;
            };

            auto seed_fallback_zero_steps = [&]() {
                fallback_step_seed_valid = true;
                fallback_s_seed.assign(n_val, 0.0);
                fallback_y_seed.assign(n_val, 0.0);
            };

            auto seed_fallback_progress = [&](const double progress_fprev_seed,
                                              const double progress_bestprog_seed,
                                              const int progress_itnfp_seed) {
                fallback_progress_seed_valid = true;
                fallback_progress_fprev = progress_fprev_seed;
                fallback_progress_bestprog = progress_bestprog_seed;
                fallback_progress_itnfp = progress_itnfp_seed;
            };

            auto seed_fallback_cg_schedule = [&](const bool valid, const double acgeps_seed,
                                                 const double bcgeps_seed, const double gpeucn20_seed,
                                                 const double gpsupn0_seed) {
                fallback_cg_schedule_seed_valid = valid;
                fallback_cg_schedule_acgeps = acgeps_seed;
                fallback_cg_schedule_bcgeps = bcgeps_seed;
                fallback_cg_schedule_gpeucn20 = gpeucn20_seed;
                fallback_cg_schedule_gpsupn0 = gpsupn0_seed;
            };

            bool tn_post_debug_captured = false;
            int tn_post_line_inform = 0;
            int tn_post_post_inform = -1;
            int tn_post_nind = 0;
            double tn_post_f_before = 0.0;
            double tn_post_f_after = 0.0;
            double tn_post_gpsupn = 0.0;
            double tn_post_gpeucn2 = 0.0;
            bool spg_post_debug_captured = false;
            int spg_post_line_inform = 0;
            int spg_post_post_inform = -1;
            int spg_post_nind = 0;
            double spg_post_f_before = 0.0;
            double spg_post_f_after = 0.0;
            double spg_post_gpsupn = 0.0;
            double spg_post_gpeucn2 = 0.0;

            std::vector<double> x_try(n_val);
            std::vector<double> g_try(n_val, 0.0);
            std::vector<int> ind_try(n_val);
            for (int i = 0; i < n_val; ++i) {
                x_try[i] = std::max(l[i], std::min(x[i], u[i]));
                ind_try[i] = i + 1;
            }

            int eval_flag = 0;
            double f_try = 0.0;
            packmol_evalal_fortran_c(n, x_try.data(), m, lambda, rho, &f_try, &eval_flag);

            if (eval_flag < 0) {
                *f = f_try;
                for (int i = 0; i < n_val; ++i) {
                    x[i] = x_try[i];
                }
                *iter = 0;
                *fcnt = 1;
                *gcnt = 0;
                *cgcnt = 0;
                *spgiter = 0;
                *spgfcnt = 0;
                *tniter = 0;
                *tnfcnt = 0;
                *tnstpcnt = 0;
                *tnintcnt = 0;
                *tnexgcnt = 0;
                *tnexbcnt = 0;
                *tnintfe = 0;
                *tnexgfe = 0;
                *tnexbfe = 0;
                *inform = eval_flag;
                return;
            }
            int precision_eval_flag = 0;
            const bool precision_solution = packmolprecision_cpp(
                n, x_try.data(), m, lambda, rho, &precision_eval_flag
            );
            if (precision_eval_flag < 0) {
                *f = f_try;
                for (int i = 0; i < n_val; ++i) {
                    x[i] = x_try[i];
                }
                *iter = 0;
                *fcnt = 1;
                *gcnt = 0;
                *cgcnt = 0;
                *spgiter = 0;
                *spgfcnt = 0;
                *tniter = 0;
                *tnfcnt = 0;
                *tnstpcnt = 0;
                *tnintcnt = 0;
                *tnexgcnt = 0;
                *tnexbcnt = 0;
                *tnintfe = 0;
                *tnexgfe = 0;
                *tnexbfe = 0;
                *inform = precision_eval_flag;
                return;
            }
            if (precision_solution) {
                *f = f_try;
                for (int i = 0; i < n_val; ++i) {
                    x[i] = x_try[i];
                }
                *iter = 0;
                *fcnt = 0;
                *gcnt = 0;
                *cgcnt = 0;
                *spgiter = 0;
                *spgfcnt = 0;
                *tniter = 0;
                *tnfcnt = 0;
                *tnstpcnt = 0;
                *tnintcnt = 0;
                *tnexgcnt = 0;
                *tnexbcnt = 0;
                *tnintfe = 0;
                *tnexgfe = 0;
                *tnexbfe = 0;
                *inform = eval_flag;
                return;
            }

            int grad_flag = 0;
            eval_gradient_full_cpp(
                n, x_try.data(), m, lambda, rho, gtype, g_try.data(), sterel, steabs, &grad_flag
            );

            if (grad_flag < 0) {
                *f = f_try;
                for (int i = 0; i < n_val; ++i) {
                    x[i] = x_try[i];
                }
                *iter = 0;
                *fcnt = 1;
                *gcnt = 1;
                *cgcnt = 0;
                *spgiter = 0;
                *spgfcnt = 0;
                *tniter = 0;
                *tnfcnt = 0;
                *tnstpcnt = 0;
                *tnintcnt = 0;
                *tnexgcnt = 0;
                *tnexbcnt = 0;
                *tnintfe = 0;
                *tnexgfe = 0;
                *tnexbfe = 0;
                *inform = grad_flag;
                return;
            }

            if (grad_flag >= 0) {
                double gpsupn_try = 0.0;
                double gpeucn2_try = 0.0;
                double gieucn2_try = 0.0;
                int nind_try = 0;
                for (int i = 0; i < n_val; ++i) {
                    const double xpg = x_try[i] - g_try[i];
                    const double gpi = std::min(u[i], std::max(l[i], xpg)) - x_try[i];
                    gpsupn_try = std::max(gpsupn_try, std::abs(gpi));
                    gpeucn2_try += gpi * gpi;
                    if (x_try[i] > l[i] && x_try[i] < u[i]) {
                        gieucn2_try += gpi * gpi;
                        ind_try[nind_try] = i + 1;
                        nind_try += 1;
                    }
                }

                const double epsgpen2 = (*epsgpen) * (*epsgpen);
                gp_ieee_signal1_cpp(
                    gpsupn_try,
                    &gencan_cg_schedule_acgeps,
                    &gencan_cg_schedule_bcgeps,
                    *cgepsf,
                    *cgepsi,
                    *cggpnf
                );
                gencan_cg_schedule_gpeucn20 = gpeucn2_try;
                gencan_cg_schedule_gpsupn0 = gpsupn_try;
                gencan_cg_schedule_seed_valid = true;
                int inform_try = -1;
                if (gpeucn2_try <= epsgpen2) {
                    inform_try = 0;
                } else if (gpsupn_try <= *epsgpsn) {
                    inform_try = 1;
                } else {
                    const double fprev = *infabs;
                    const double currprog = fprev - f_try;
                    const double bestprog = std::max(currprog, 0.0);
                    int itnfp = 0;
                    if (currprog <= (*epsnfp) * bestprog) {
                        itnfp += 1;
                        if (itnfp >= *maxitnfp) {
                            inform_try = 2;
                        }
                    }
                }

                if (inform_try < 0) {
                    double gpnmax = *infabs;
                    if (gpeucn2_try >= gpnmax) {
                        inform_try = 3;
                    }
                }

                if (inform_try < 0 && f_try <= *fmin) {
                    inform_try = 4;
                } else if (inform_try < 0 && 0 >= *maxit) {
                    inform_try = 7;
                } else if (inform_try < 0 && 1 >= *maxfc) {
                    inform_try = 8;
                }

                    if (inform_try >= 0) {
                        *f = f_try;
                    for (int i = 0; i < n_val; ++i) {
                        x[i] = x_try[i];
                        g[i] = g_try[i];
                    }
                    for (int i = 0; i < nind_try; ++i) {
                        ind[i] = ind_try[i];
                    }
                    *gpeucn2 = gpeucn2_try;
                    *gpsupn = gpsupn_try;
                    *iter = 0;
                    *fcnt = 1;
                    *gcnt = 1;
                    *cgcnt = 0;
                    *spgiter = 0;
                    *spgfcnt = 0;
                    *tniter = 0;
                    *tnfcnt = 0;
                    *tnstpcnt = 0;
                    *tnintcnt = 0;
                    *tnexgcnt = 0;
                    *tnexbcnt = 0;
                        *tnintfe = 0;
                        *tnexgfe = 0;
                        *tnexbfe = 0;
                        if (*maxitngp > 0 && inform_try == 3) {
                            lastgpns[0] = gpeucn2_try;
                        }
                        *inform = inform_try;
                        return;
                    }

                double progress_fprev = f_try;
                double progress_bestprog = std::max((*infabs) - f_try, 0.0);
                int progress_itnfp = 0;

                const double ometa2 = (1.0 - *eta) * (1.0 - *eta);
                if (gpeucn2_try > 0.0 && gieucn2_try <= ometa2 * gpeucn2_try) {
                    std::vector<double> x_work = x_try;
                    std::vector<double> g_work = g_try;
                    std::vector<double> x_prev = x_try;
                    std::vector<double> g_prev = g_try;
                    std::vector<double> s_work(n_val, 0.0);
                    std::vector<double> y_work(n_val, 0.0);
                    std::vector<double> xtrial_work(n_val, 0.0);
                    std::vector<double> d_work(n_val, 0.0);
                    std::vector<int> ind_work(n_val, 0);

                    int iter_work = 1;
                    int fcnt_work = 1;
                    int gcnt_work = 1;
                    int cgcnt_work = 0;
                    int spgiter_work = 1;
                    int spgfcnt_work = 0;
                    int tniter_work = 0;
                    int tnfcnt_work = 0;
                    int tnstpcnt_work = 0;
                    int tnintcnt_work = 0;
                    int tnexgcnt_work = 0;
                    int tnexbcnt_work = 0;
                    int tnintfe_work = 0;
                    int tnexgfe_work = 0;
                    int tnexbfe_work = 0;

                    double f_work = f_try;

                    double xnorm = 0.0;
                    for (int i = 0; i < n_val; ++i) {
                        xnorm += x_work[i] * x_work[i];
                    }
                    xnorm = std::sqrt(xnorm);

                    double lamspg = std::max(1.0, xnorm) / std::sqrt(gpeucn2_try);
                    lamspg = std::min(*lspgma, std::max(*lspgmi, lamspg));

                    int ls_inform = 0;
                    const int fcnt_prev = fcnt_work;
                    spgls_cpp(
                        n, x_work.data(), m, lambda, rho, &f_work, g_work.data(), l, u, &lamspg,
                        nint, mininterp, fmin, maxfc, &fcnt_work, &ls_inform, xtrial_work.data(),
                        d_work.data(), gamma, sigma1, sigma2, epsrel, epsabs
                    );
                    spgfcnt_work += (fcnt_work - fcnt_prev);

                    if (ls_inform < 0) {
                        *f = f_work;
                        for (int i = 0; i < n_val; ++i) {
                            x[i] = x_work[i];
                            g[i] = g_work[i];
                        }
                        *iter = iter_work;
                        *fcnt = fcnt_work;
                        *gcnt = gcnt_work;
                        *cgcnt = cgcnt_work;
                        *spgiter = spgiter_work;
                        *spgfcnt = spgfcnt_work;
                        *tniter = tniter_work;
                        *tnfcnt = tnfcnt_work;
                        *tnstpcnt = tnstpcnt_work;
                        *tnintcnt = tnintcnt_work;
                        *tnexgcnt = tnexgcnt_work;
                        *tnexbcnt = tnexbcnt_work;
                        *tnintfe = tnintfe_work;
                        *tnexgfe = tnexgfe_work;
                        *tnexbfe = tnexbfe_work;
                        *inform = ls_inform;
                        return;
                    }

                    int grad_after_spg = 0;
                    eval_gradient_full_cpp(
                        n, x_work.data(), m, lambda, rho, gtype, g_work.data(), sterel, steabs, &grad_after_spg
                    );
                    gcnt_work += 1;

                    if (grad_after_spg < 0) {
                        *f = f_work;
                        for (int i = 0; i < n_val; ++i) {
                            x[i] = x_work[i];
                            g[i] = g_work[i];
                        }
                        *iter = iter_work;
                        *fcnt = fcnt_work;
                        *gcnt = gcnt_work;
                        *cgcnt = cgcnt_work;
                        *spgiter = spgiter_work;
                        *spgfcnt = spgfcnt_work;
                        *tniter = tniter_work;
                        *tnfcnt = tnfcnt_work;
                        *tnstpcnt = tnstpcnt_work;
                        *tnintcnt = tnintcnt_work;
                        *tnexgcnt = tnexgcnt_work;
                        *tnexbcnt = tnexbcnt_work;
                        *tnintfe = tnintfe_work;
                        *tnexgfe = tnexgfe_work;
                        *tnexbfe = tnexbfe_work;
                        *inform = grad_after_spg;
                        return;
                    }

                    for (int i = 0; i < n_val; ++i) {
                        if (x_work[i] <= l[i] + std::max((*epsrel) * std::abs(l[i]), *epsabs)) {
                            x_work[i] = l[i];
                        } else if (x_work[i] >= u[i] - std::max((*epsrel) * std::abs(u[i]), *epsabs)) {
                            x_work[i] = u[i];
                        }
                    }

                    double gpsupn_after = 0.0;
                    double gpeucn2_after = 0.0;
                    int nind_after = 0;
                    for (int i = 0; i < n_val; ++i) {
                        s_work[i] = x_work[i] - x_prev[i];
                        y_work[i] = g_work[i] - g_prev[i];
                        const double xpg = x_work[i] - g_work[i];
                        const double gpi = std::min(u[i], std::max(l[i], xpg)) - x_work[i];
                        gpsupn_after = std::max(gpsupn_after, std::abs(gpi));
                        gpeucn2_after += gpi * gpi;
                        if (x_work[i] > l[i] && x_work[i] < u[i]) {
                            ind_work[nind_after] = i + 1;
                            nind_after += 1;
                        }
                    }

                    if (ls_inform == 6) {
                        *f = f_work;
                        for (int i = 0; i < n_val; ++i) {
                            x[i] = x_work[i];
                            g[i] = g_work[i];
                            s[i] = s_work[i];
                            y[i] = y_work[i];
                            d[i] = d_work[i];
                        }
                        for (int i = 0; i < nind_after; ++i) {
                            ind[i] = ind_work[i];
                        }
                        *gpeucn2 = gpeucn2_after;
                        *gpsupn = gpsupn_after;
                        *iter = iter_work;
                        *fcnt = fcnt_work;
                        *gcnt = gcnt_work;
                        *cgcnt = cgcnt_work;
                        *spgiter = spgiter_work;
                        *spgfcnt = spgfcnt_work;
                        *tniter = tniter_work;
                        *tnfcnt = tnfcnt_work;
                        *tnstpcnt = tnstpcnt_work;
                        *tnintcnt = tnintcnt_work;
                        *tnexgcnt = tnexgcnt_work;
                        *tnexbcnt = tnexbcnt_work;
                        *tnintfe = tnintfe_work;
                        *tnexgfe = tnexgfe_work;
                        *tnexbfe = tnexbfe_work;
                        *inform = 6;
                        return;
                    }

                    int precision_after_spg_flag = 0;
                    const bool precision_after_spg = packmolprecision_cpp(
                        n, x_work.data(), m, lambda, rho, &precision_after_spg_flag
                    );
                    if (precision_after_spg_flag < 0) {
                        *f = f_work;
                        for (int i = 0; i < n_val; ++i) {
                            x[i] = x_work[i];
                            g[i] = g_work[i];
                        }
                        *iter = iter_work;
                        *fcnt = fcnt_work;
                        *gcnt = gcnt_work;
                        *cgcnt = cgcnt_work;
                        *spgiter = spgiter_work;
                        *spgfcnt = spgfcnt_work;
                        *tniter = tniter_work;
                        *tnfcnt = tnfcnt_work;
                        *tnstpcnt = tnstpcnt_work;
                        *tnintcnt = tnintcnt_work;
                        *tnexgcnt = tnexgcnt_work;
                        *tnexbcnt = tnexbcnt_work;
                        *tnintfe = tnintfe_work;
                        *tnexgfe = tnexgfe_work;
                        *tnexbfe = tnexbfe_work;
                        *inform = precision_after_spg_flag;
                        return;
                    }
                    int post_inform = evaluate_post_step_inform_cpp(
                        precision_after_spg,
                        ls_inform,
                        gpeucn2_after,
                        *epsgpen,
                        gpsupn_after,
                        *epsgpsn,
                        f_try,
                        f_work,
                        *epsnfp,
                        *maxitnfp,
                        *infabs,
                        lastgpns,
                        *maxitngp,
                        *fmin,
                        iter_work,
                        *maxit,
                        fcnt_work,
                        *maxfc,
                        &progress_fprev,
                        &progress_bestprog,
                        &progress_itnfp
                    );

                    int retry_budget = spg_post_retry_steps();
                    while (post_inform < 0 && retry_budget > 0) {
                        retry_budget -= 1;
                        spgiter_work += 1;
                        const double xnorm_for_spg = std::sqrt(norm2_kernel(n_val, x_work.data()));
                        double lamspg_retry = std::max(1.0, xnorm_for_spg) / std::sqrt(std::max(gpeucn2_after, 1.0e-30));
                        lamspg_retry = std::min(*lspgma, std::max(*lspgmi, lamspg_retry));

                        const double f_before_retry = f_work;
                        int spg_line_inform = 0;
                        const int fcnt_prev_retry = fcnt_work;
                        spgls_cpp(
                            n, x_work.data(), m, lambda, rho, &f_work, g_work.data(), l, u, &lamspg_retry,
                            nint, mininterp, fmin, maxfc, &fcnt_work, &spg_line_inform, xtrial_work.data(),
                            d_work.data(), gamma, sigma1, sigma2, epsrel, epsabs
                        );
                        spgfcnt_work += (fcnt_work - fcnt_prev_retry);

                        if (spg_line_inform < 0) {
                            *f = f_work;
                            for (int i = 0; i < n_val; ++i) {
                                x[i] = x_work[i];
                                g[i] = g_work[i];
                            }
                            *iter = iter_work;
                            *fcnt = fcnt_work;
                            *gcnt = gcnt_work;
                            *cgcnt = cgcnt_work;
                            *spgiter = spgiter_work;
                            *spgfcnt = spgfcnt_work;
                            *tniter = tniter_work;
                            *tnfcnt = tnfcnt_work;
                            *tnstpcnt = tnstpcnt_work;
                            *tnintcnt = tnintcnt_work;
                            *tnexgcnt = tnexgcnt_work;
                            *tnexbcnt = tnexbcnt_work;
                            *tnintfe = tnintfe_work;
                            *tnexgfe = tnexgfe_work;
                            *tnexbfe = tnexbfe_work;
                            *inform = spg_line_inform;
                            return;
                        }

                        int grad_after_retry = 0;
                        eval_gradient_full_cpp(
                            n, x_work.data(), m, lambda, rho, gtype, g_work.data(), sterel, steabs, &grad_after_retry
                        );
                        gcnt_work += 1;

                        if (grad_after_retry < 0) {
                            *f = f_work;
                            for (int i = 0; i < n_val; ++i) {
                                x[i] = x_work[i];
                                g[i] = g_work[i];
                            }
                            *iter = iter_work;
                            *fcnt = fcnt_work;
                            *gcnt = gcnt_work;
                            *cgcnt = cgcnt_work;
                            *spgiter = spgiter_work;
                            *spgfcnt = spgfcnt_work;
                            *tniter = tniter_work;
                            *tnfcnt = tnfcnt_work;
                            *tnstpcnt = tnstpcnt_work;
                            *tnintcnt = tnintcnt_work;
                            *tnexgcnt = tnexgcnt_work;
                            *tnexbcnt = tnexbcnt_work;
                            *tnintfe = tnintfe_work;
                            *tnexgfe = tnexgfe_work;
                            *tnexbfe = tnexbfe_work;
                            *inform = grad_after_retry;
                            return;
                        }

                        gpsupn_after = 0.0;
                        gpeucn2_after = 0.0;
                        nind_after = 0;
                        for (int i = 0; i < n_val; ++i) {
                            if (x_work[i] <= l[i] + std::max((*epsrel) * std::abs(l[i]), *epsabs)) {
                                x_work[i] = l[i];
                            } else if (x_work[i] >= u[i] - std::max((*epsrel) * std::abs(u[i]), *epsabs)) {
                                x_work[i] = u[i];
                            }
                            s_work[i] = x_work[i] - x_prev[i];
                            y_work[i] = g_work[i] - g_prev[i];
                            const double xpg = x_work[i] - g_work[i];
                            const double gpi = std::min(u[i], std::max(l[i], xpg)) - x_work[i];
                            gpsupn_after = std::max(gpsupn_after, std::abs(gpi));
                            gpeucn2_after += gpi * gpi;
                            if (x_work[i] > l[i] && x_work[i] < u[i]) {
                                ind_work[nind_after] = i + 1;
                                nind_after += 1;
                            }
                        }

                        int precision_after_retry_flag = 0;
                        const bool precision_after_retry = packmolprecision_cpp(
                            n, x_work.data(), m, lambda, rho, &precision_after_retry_flag
                        );
                        if (precision_after_retry_flag < 0) {
                            *f = f_work;
                            for (int i = 0; i < n_val; ++i) {
                                x[i] = x_work[i];
                                g[i] = g_work[i];
                            }
                            *iter = iter_work;
                            *fcnt = fcnt_work;
                            *gcnt = gcnt_work;
                            *cgcnt = cgcnt_work;
                            *spgiter = spgiter_work;
                            *spgfcnt = spgfcnt_work;
                            *tniter = tniter_work;
                            *tnfcnt = tnfcnt_work;
                            *tnstpcnt = tnstpcnt_work;
                            *tnintcnt = tnintcnt_work;
                            *tnexgcnt = tnexgcnt_work;
                            *tnexbcnt = tnexbcnt_work;
                            *tnintfe = tnintfe_work;
                            *tnexgfe = tnexgfe_work;
                            *tnexbfe = tnexbfe_work;
                            *inform = precision_after_retry_flag;
                            return;
                        }
                        ls_inform = spg_line_inform;
                        post_inform = evaluate_post_step_inform_cpp(
                            precision_after_retry,
                            ls_inform,
                            gpeucn2_after,
                            *epsgpen,
                            gpsupn_after,
                            *epsgpsn,
                            f_before_retry,
                            f_work,
                            *epsnfp,
                            *maxitnfp,
                            *infabs,
                            lastgpns,
                            *maxitngp,
                            *fmin,
                            iter_work,
                            *maxit,
                            fcnt_work,
                            *maxfc,
                            &progress_fprev,
                            &progress_bestprog,
                            &progress_itnfp
                        );

                        if (gencan_debug_enabled()) {
                            std::fprintf(
                                stderr,
                                "[gencan-cpp-spg-retry] step_left=%d line_inform=%d post_inform=%d f=%.16e gpsupn=%.16e gpeucn2=%.16e\n",
                                retry_budget,
                                ls_inform,
                                post_inform,
                                f_work,
                                gpsupn_after,
                                gpeucn2_after
                            );
                        }
                    }

                    spg_post_debug_captured = true;
                    spg_post_line_inform = ls_inform;
                    spg_post_post_inform = post_inform;
                    spg_post_nind = nind_after;
                    spg_post_f_before = f_try;
                    spg_post_f_after = f_work;
                    spg_post_gpsupn = gpsupn_after;
                    spg_post_gpeucn2 = gpeucn2_after;

                    if (post_inform >= 0) {
                        *f = f_work;
                        for (int i = 0; i < n_val; ++i) {
                            x[i] = x_work[i];
                            g[i] = g_work[i];
                            s[i] = s_work[i];
                            y[i] = y_work[i];
                            d[i] = d_work[i];
                        }
                        for (int i = 0; i < nind_after; ++i) {
                            ind[i] = ind_work[i];
                        }
                        *gpeucn2 = gpeucn2_after;
                        *gpsupn = gpsupn_after;
                        *iter = iter_work;
                        *fcnt = fcnt_work;
                        *gcnt = gcnt_work;
                        *cgcnt = cgcnt_work;
                        *spgiter = spgiter_work;
                        *spgfcnt = spgfcnt_work;
                        *tniter = tniter_work;
                        *tnfcnt = tnfcnt_work;
                        *tnstpcnt = tnstpcnt_work;
                        *tnintcnt = tnintcnt_work;
                        *tnexgcnt = tnexgcnt_work;
                        *tnexbcnt = tnexbcnt_work;
                        *tnintfe = tnintfe_work;
                        *tnexgfe = tnexgfe_work;
                        *tnexbfe = tnexbfe_work;
                        if (*maxitngp > 0) {
                            lastgpns[iter_work % (*maxitngp)] = gpeucn2_after;
                        }
                        *inform = post_inform;
                        return;
                    }
                    if (!spg_post_handoff_enabled()) {
                        *f = f_work;
                        for (int i = 0; i < n_val; ++i) {
                            x[i] = x_work[i];
                            g[i] = g_work[i];
                            s[i] = s_work[i];
                            y[i] = y_work[i];
                            d[i] = d_work[i];
                        }
                        for (int i = 0; i < nind_after; ++i) {
                            ind[i] = ind_work[i];
                        }
                        *gpeucn2 = gpeucn2_after;
                        *gpsupn = gpsupn_after;
                        *iter = iter_work;
                        *fcnt = fcnt_work;
                        *gcnt = gcnt_work;
                        *cgcnt = cgcnt_work;
                        *spgiter = spgiter_work;
                        *spgfcnt = spgfcnt_work;
                        *tniter = tniter_work;
                        *tnfcnt = tnfcnt_work;
                        *tnstpcnt = tnstpcnt_work;
                        *tnintcnt = tnintcnt_work;
                        *tnexgcnt = tnexgcnt_work;
                        *tnexbcnt = tnexbcnt_work;
                        *tnintfe = tnintfe_work;
                        *tnexgfe = tnexgfe_work;
                        *tnexbfe = tnexbfe_work;
                        if (*maxitngp > 0) {
                            lastgpns[iter_work % (*maxitngp)] = gpeucn2_after;
                        }

                        if (gencan_debug_enabled()) {
                            std::fprintf(
                                stderr,
                                "[gencan-cpp-spg-post-direct-tail] mode=%d iter=%d fcnt=%d gcnt=%d\n",
                                static_cast<int>(active_impl_mode()),
                                *iter,
                                *fcnt,
                                *gcnt
                            );
                        }

                        continue_spg_post_nonterminal_cpp(
                            n_val,
                            n,
                            x,
                            l,
                            u,
                            m,
                            lambda,
                            rho,
                            epsgpen,
                            epsgpsn,
                            maxitnfp,
                            epsnfp,
                            maxitngp,
                            fmin,
                            maxit,
                            maxfc,
                            udelta0,
                            ucgmaxit,
                            cgscre,
                            cggpnf,
                            cgepsi,
                            cgepsf,
                            epsnqmp,
                            maxitnqmp,
                            nearlyq,
                            nint,
                            next,
                            mininterp,
                            maxextrap,
                            gtype,
                            htvtype,
                            trtype,
                            iprint,
                            ncomp,
                            f,
                            g,
                            gpeucn2,
                            gpsupn,
                            iter,
                            fcnt,
                            gcnt,
                            cgcnt,
                            spgiter,
                            spgfcnt,
                            tniter,
                            tnfcnt,
                            tnstpcnt,
                            tnintcnt,
                            tnexgcnt,
                            tnexbcnt,
                            tnintfe,
                            tnexgfe,
                            tnexbfe,
                            inform,
                            s,
                            y,
                            d,
                            ind,
                            lastgpns,
                            w,
                            eta,
                            delmin,
                            lspgma,
                            lspgmi,
                            theta,
                            gamma,
                            beta,
                            sigma1,
                            sigma2,
                            sterel,
                            steabs,
                            epsrel,
                            epsabs,
                            infrel,
                            infabs,
                            progress_fprev,
                            progress_bestprog,
                            progress_itnfp
                        );
                        return;
                    }

                    fallback_reason = "spg_post_nonterminal";
                    seed_fallback_state("spg_post_nonterminal", x_work, &g_work, f_work);
                    seed_fallback_counters(
                        iter_work, fcnt_work, gcnt_work, cgcnt_work, spgiter_work, spgfcnt_work,
                        tniter_work, tnfcnt_work, tnstpcnt_work, tnintcnt_work, tnexgcnt_work,
                        tnexbcnt_work, tnintfe_work, tnexgfe_work, tnexbfe_work
                    );
                    seed_fallback_steps(s_work, y_work);
                    seed_fallback_progress(progress_fprev, progress_bestprog, progress_itnfp);
                } else {
                    std::vector<double> x_work = x_try;
                    std::vector<double> g_work = g_try;
                    std::vector<double> l_work(n_val);
                    std::vector<double> u_work(n_val);
                    for (int i = 0; i < n_val; ++i) {
                        l_work[i] = l[i];
                        u_work[i] = u[i];
                    }

                    const int nind_work = nind_try;
                    if (nind_work > 0) {
                        shrink_inplace(nind_work, ind_try.data(), x_work.data());
                        shrink_inplace(nind_work, ind_try.data(), g_work.data());
                        shrink_inplace(nind_work, ind_try.data(), l_work.data());
                        shrink_inplace(nind_work, ind_try.data(), u_work.data());

                        int iter_work = 1;
                        int fcnt_work = 1;
                        int gcnt_work = 1;
                        int cgcnt_work = 0;
                        int spgiter_work = 0;
                        int spgfcnt_work = 0;
                        int tniter_work = 1;
                        int tnfcnt_work = 0;
                        int tnstpcnt_work = 0;
                        int tnintcnt_work = 0;
                        int tnexgcnt_work = 0;
                        int tnexbcnt_work = 0;
                        int tnintfe_work = 0;
                        int tnexgfe_work = 0;
                        int tnexbfe_work = 0;
                        const int pre_tn_iter_seed = iter_work;
                        const int pre_tn_fcnt_seed = fcnt_work;
                        const int pre_tn_gcnt_seed = gcnt_work;
                        const int pre_tn_cgcnt_seed = cgcnt_work;
                        const int pre_tn_spgiter_seed = spgiter_work;
                        const int pre_tn_spgfcnt_seed = spgfcnt_work;
                        const int pre_tn_tniter_seed = tniter_work;
                        const int pre_tn_tnfcnt_seed = tnfcnt_work;
                        const int pre_tn_tnstpcnt_seed = tnstpcnt_work;
                        const int pre_tn_tnintcnt_seed = tnintcnt_work;
                        const int pre_tn_tnexgcnt_seed = tnexgcnt_work;
                        const int pre_tn_tnexbcnt_seed = tnexbcnt_work;
                        const int pre_tn_tnintfe_seed = tnintfe_work;
                        const int pre_tn_tnexgfe_seed = tnexgfe_work;
                        const int pre_tn_tnexbfe_seed = tnexbfe_work;

                        const double xnorm_try = std::sqrt(norm2_kernel(n_val, x_try.data()));
                        double delta = 0.0;
                        if (*udelta0 <= 0.0) {
                            delta = std::max(*delmin, 0.1 * std::max(1.0, xnorm_try));
                        } else {
                            delta = *udelta0;
                        }

                        const double epsgpen2 = (*epsgpen) * (*epsgpen);
                        double acgeps = 0.0;
                        double bcgeps = 0.0;
                        if (gencan_cg_schedule_seed_valid) {
                            acgeps = gencan_cg_schedule_acgeps;
                            bcgeps = gencan_cg_schedule_bcgeps;
                        } else {
                            gp_ieee_signal1_cpp(
                                gpsupn_try, &acgeps, &bcgeps, *cgepsf, *cgepsi, *cggpnf
                            );
                        }
                        const double gpeucn20 = gencan_cg_schedule_seed_valid
                            ? gencan_cg_schedule_gpeucn20
                            : gpeucn2_try;
                        const double gpsupn0 = gencan_cg_schedule_seed_valid
                            ? gencan_cg_schedule_gpsupn0
                            : gpsupn_try;
                        double kappa = 0.0;
                        double cgeps = *cgepsf;
                        int cgmaxit = 0;
                        gp_ieee_signal2_cpp(
                            &cgmaxit, nind_work, *nearlyq, *ucgmaxit, *cgscre, &kappa, gpeucn2_try,
                            gpeucn20, epsgpen2, *epsgpsn, &cgeps, acgeps, bcgeps, *cgepsf, *cgepsi,
                            gpsupn_try, gpsupn0
                        );

                        std::vector<double> cg_s(n_val, 0.0);
                        std::vector<double> cg_w(n_val, 0.0);
                        std::vector<double> cg_y(n_val, 0.0);
                        std::vector<double> cg_r(n_val, 0.0);
                        std::vector<double> cg_d(n_val, 0.0);
                        std::vector<double> cg_sprev(n_val, 0.0);

                        int cg_iter = 0;
                        int rbdtype = 0;
                        int rbdind = 0;
                        int cg_inform = 0;
                        packmol_gencan_cg_bridge(
                            &nind_work, ind_try.data(), n, x_work.data(), m, lambda, rho, g_work.data(),
                            &delta, l_work.data(), u_work.data(), &cgeps, epsnqmp, maxitnqmp, &cgmaxit,
                            nearlyq, gtype, htvtype, trtype, iprint, ncomp, cg_s.data(), &cg_iter, &rbdtype,
                            &rbdind, &cg_inform, cg_w.data(), cg_y.data(), cg_r.data(), cg_d.data(),
                            cg_sprev.data(), theta, sterel, steabs, epsrel, epsabs, infrel, infabs
                        );
                        cgcnt_work += cg_iter;

                        if (cg_inform < 0) {
                            *f = f_try;
                            for (int i = 0; i < n_val; ++i) {
                                x[i] = x_work[i];
                                g[i] = g_work[i];
                            }
                            *iter = iter_work;
                            *fcnt = fcnt_work;
                            *gcnt = gcnt_work;
                            *cgcnt = cgcnt_work;
                            *spgiter = spgiter_work;
                            *spgfcnt = spgfcnt_work;
                            *tniter = tniter_work;
                            *tnfcnt = tnfcnt_work;
                            *tnstpcnt = tnstpcnt_work;
                            *tnintcnt = tnintcnt_work;
                            *tnexgcnt = tnexgcnt_work;
                            *tnexbcnt = tnexbcnt_work;
                            *tnintfe = tnintfe_work;
                            *tnexgfe = tnexgfe_work;
                            *tnexbfe = tnexbfe_work;
                            *inform = cg_inform;
                            return;
                        }

                        if (cg_inform >= 0) {
                            double amax = *infabs;
                            if (cg_inform == 2) {
                                amax = 1.0;
                            } else {
                                for (int i = 0; i < nind_work; ++i) {
                                    if (cg_s[i] > 0.0) {
                                        const double amaxx = (u_work[i] - x_work[i]) / cg_s[i];
                                        if (amaxx < amax) {
                                            amax = amaxx;
                                            rbdind = i + 1;
                                            rbdtype = 2;
                                        }
                                    } else if (cg_s[i] < 0.0) {
                                        const double amaxx = (l_work[i] - x_work[i]) / cg_s[i];
                                        if (amaxx < amax) {
                                            amax = amaxx;
                                            rbdind = i + 1;
                                            rbdtype = 1;
                                        }
                                    }
                                }
                            }

                            int tnls_inform = 0;
                            const int tnint_prev = tnintcnt_work;
                            const int tnexg_prev = tnexgcnt_work;
                            const int tnexb_prev = tnexbcnt_work;
                            const int fcnt_prev = fcnt_work;
                            double f_work = f_try;
                            std::vector<double> xplus(n_val, 0.0);
                            std::vector<double> xtmp(n_val, 0.0);
                            std::vector<double> xbext(n_val, 0.0);
                            packmol_gencan_tnls_bridge(
                                &nind_work, ind_try.data(), n, x_work.data(), m, lambda, rho, l_work.data(),
                                u_work.data(), &f_work, g_work.data(), cg_s.data(), &amax, &rbdtype, &rbdind,
                                nint, next, mininterp, maxextrap, fmin, maxfc, gtype, iprint, &fcnt_work,
                                &gcnt_work, &tnintcnt_work, &tnexgcnt_work, &tnexbcnt_work, &tnls_inform,
                                xplus.data(), xtmp.data(), xbext.data(), gamma, beta, sigma1, sigma2, sterel,
                                steabs, epsrel, epsabs, infrel, infabs
                            );

                            if (tnls_inform < 0) {
                                *f = f_work;
                                for (int i = 0; i < n_val; ++i) {
                                    x[i] = x_work[i];
                                    g[i] = g_work[i];
                                }
                                *iter = iter_work;
                                *fcnt = fcnt_work;
                                *gcnt = gcnt_work;
                                *cgcnt = cgcnt_work;
                                *spgiter = spgiter_work;
                                *spgfcnt = spgfcnt_work;
                                *tniter = tniter_work;
                                *tnfcnt = tnfcnt_work;
                                *tnstpcnt = tnstpcnt_work;
                                *tnintcnt = tnintcnt_work;
                                *tnexgcnt = tnexgcnt_work;
                                *tnexbcnt = tnexbcnt_work;
                                *tnintfe = tnintfe_work;
                                *tnexgfe = tnexgfe_work;
                                *tnexbfe = tnexbfe_work;
                                *inform = tnls_inform;
                                return;
                            }

                            if (tnls_inform >= 0) {
                                if (tnintcnt_work > tnint_prev) {
                                    tnintfe_work += (fcnt_work - fcnt_prev);
                                } else if (tnexgcnt_work > tnexg_prev) {
                                    tnexgfe_work += (fcnt_work - fcnt_prev);
                                } else if (tnexbcnt_work > tnexb_prev) {
                                    tnexbfe_work += (fcnt_work - fcnt_prev);
                                } else {
                                    tnstpcnt_work += 1;
                                }
                                tnfcnt_work += (fcnt_work - fcnt_prev);

                                expand_inplace(nind_work, ind_try.data(), x_work.data());
                                expand_inplace(nind_work, ind_try.data(), g_work.data());
                                expand_inplace(nind_work, ind_try.data(), l_work.data());
                                expand_inplace(nind_work, ind_try.data(), u_work.data());

                                for (int i = 0; i < n_val; ++i) {
                                    if (x_work[i] <= l_work[i] + std::max((*epsrel) * std::abs(l_work[i]), *epsabs)) {
                                        x_work[i] = l_work[i];
                                    } else if (x_work[i] >= u_work[i] - std::max((*epsrel) * std::abs(u_work[i]), *epsabs)) {
                                        x_work[i] = u_work[i];
                                    }
                                }

                                int line_inform = tnls_inform;
                                if (tnls_inform == 6) {
                                    spgiter_work += 1;
                                    double xnorm_for_spg = std::sqrt(norm2_kernel(n_val, x_work.data()));
                                    double lamspg = std::max(1.0, xnorm_for_spg) / std::sqrt(gpeucn2_try);
                                    lamspg = std::min(*lspgma, std::max(*lspgmi, lamspg));

                                    const int fcnt_prev_spg = fcnt_work;
                                    spgls_cpp(
                                        n, x_work.data(), m, lambda, rho, &f_work, g_work.data(), l_work.data(),
                                        u_work.data(), &lamspg, nint, mininterp, fmin, maxfc, &fcnt_work,
                                        &line_inform, xplus.data(), xtmp.data(), gamma, sigma1, sigma2, epsrel,
                                        epsabs
                                    );
                                    spgfcnt_work += (fcnt_work - fcnt_prev_spg);

                                    if (line_inform < 0) {
                                        *f = f_work;
                                        for (int i = 0; i < n_val; ++i) {
                                            x[i] = x_work[i];
                                            g[i] = g_work[i];
                                        }
                                        *iter = iter_work;
                                        *fcnt = fcnt_work;
                                        *gcnt = gcnt_work;
                                        *cgcnt = cgcnt_work;
                                        *spgiter = spgiter_work;
                                        *spgfcnt = spgfcnt_work;
                                        *tniter = tniter_work;
                                        *tnfcnt = tnfcnt_work;
                                        *tnstpcnt = tnstpcnt_work;
                                        *tnintcnt = tnintcnt_work;
                                        *tnexgcnt = tnexgcnt_work;
                                        *tnexbcnt = tnexbcnt_work;
                                        *tnintfe = tnintfe_work;
                                        *tnexgfe = tnexgfe_work;
                                        *tnexbfe = tnexbfe_work;
                                        *inform = line_inform;
                                        return;
                                    }

                                    int grad_after_spg = 0;
                                    eval_gradient_full_cpp(
                                        n, x_work.data(), m, lambda, rho, gtype, g_work.data(), sterel, steabs, &grad_after_spg
                                    );
                                    gcnt_work += 1;

                                    if (grad_after_spg < 0) {
                                        *f = f_work;
                                        for (int i = 0; i < n_val; ++i) {
                                            x[i] = x_work[i];
                                            g[i] = g_work[i];
                                        }
                                        *iter = iter_work;
                                        *fcnt = fcnt_work;
                                        *gcnt = gcnt_work;
                                        *cgcnt = cgcnt_work;
                                        *spgiter = spgiter_work;
                                        *spgfcnt = spgfcnt_work;
                                        *tniter = tniter_work;
                                        *tnfcnt = tnfcnt_work;
                                        *tnstpcnt = tnstpcnt_work;
                                        *tnintcnt = tnintcnt_work;
                                        *tnexgcnt = tnexgcnt_work;
                                        *tnexbcnt = tnexbcnt_work;
                                        *tnintfe = tnintfe_work;
                                        *tnexgfe = tnexgfe_work;
                                        *tnexbfe = tnexbfe_work;
                                        *inform = grad_after_spg;
                                        return;
                                    }
                                }

                                double gpsupn_after = 0.0;
                                double gpeucn2_after = 0.0;
                                for (int i = 0; i < n_val; ++i) {
                                    const double xpg = x_work[i] - g_work[i];
                                    const double gpi = std::min(u_work[i], std::max(l_work[i], xpg)) - x_work[i];
                                    gpsupn_after = std::max(gpsupn_after, std::abs(gpi));
                                    gpeucn2_after += gpi * gpi;
                                }

                                int precision_after_flag = 0;
                                const bool precision_after = packmolprecision_cpp(
                                    n, x_work.data(), m, lambda, rho, &precision_after_flag
                                );
                                if (precision_after_flag < 0) {
                                    *f = f_work;
                                    for (int i = 0; i < n_val; ++i) {
                                        x[i] = x_work[i];
                                        g[i] = g_work[i];
                                    }
                                    *iter = iter_work;
                                    *fcnt = fcnt_work;
                                    *gcnt = gcnt_work;
                                    *cgcnt = cgcnt_work;
                                    *spgiter = spgiter_work;
                                    *spgfcnt = spgfcnt_work;
                                    *tniter = tniter_work;
                                    *tnfcnt = tnfcnt_work;
                                    *tnstpcnt = tnstpcnt_work;
                                    *tnintcnt = tnintcnt_work;
                                    *tnexgcnt = tnexgcnt_work;
                                    *tnexbcnt = tnexbcnt_work;
                                    *tnintfe = tnintfe_work;
                                    *tnexgfe = tnexgfe_work;
                                    *tnexbfe = tnexbfe_work;
                                    *inform = precision_after_flag;
                                    return;
                                }
                                const double pre_tn_progress_fprev = progress_fprev;
                                const double pre_tn_progress_bestprog = progress_bestprog;
                                const int pre_tn_progress_itnfp = progress_itnfp;
                                int post_inform = evaluate_post_step_inform_cpp(
                                    precision_after,
                                    line_inform,
                                    gpeucn2_after,
                                    *epsgpen,
                                    gpsupn_after,
                                    *epsgpsn,
                                    f_try,
                                    f_work,
                                    *epsnfp,
                                    *maxitnfp,
                                    *infabs,
                                    lastgpns,
                                    *maxitngp,
                                    *fmin,
                                    iter_work,
                                    *maxit,
                                    fcnt_work,
                                    *maxfc,
                                    &progress_fprev,
                                    &progress_bestprog,
                                    &progress_itnfp
                                );

                                int retry_budget = tn_post_retry_spg_steps();
                                while (post_inform < 0 && retry_budget > 0) {
                                    retry_budget -= 1;
                                    spgiter_work += 1;
                                    const double xnorm_for_spg = std::sqrt(norm2_kernel(n_val, x_work.data()));
                                    double lamspg_retry = std::max(1.0, xnorm_for_spg) / std::sqrt(std::max(gpeucn2_after, 1.0e-30));
                                    lamspg_retry = std::min(*lspgma, std::max(*lspgmi, lamspg_retry));

                                    const double f_before_retry = f_work;
                                    int spg_line_inform = 0;
                                    const int fcnt_prev_spg = fcnt_work;
                                    spgls_cpp(
                                        n, x_work.data(), m, lambda, rho, &f_work, g_work.data(), l_work.data(),
                                        u_work.data(), &lamspg_retry, nint, mininterp, fmin, maxfc, &fcnt_work,
                                        &spg_line_inform, xplus.data(), xtmp.data(), gamma, sigma1, sigma2, epsrel,
                                        epsabs
                                    );
                                    spgfcnt_work += (fcnt_work - fcnt_prev_spg);

                                    if (spg_line_inform < 0) {
                                        *f = f_work;
                                        for (int i = 0; i < n_val; ++i) {
                                            x[i] = x_work[i];
                                            g[i] = g_work[i];
                                        }
                                        *iter = iter_work;
                                        *fcnt = fcnt_work;
                                        *gcnt = gcnt_work;
                                        *cgcnt = cgcnt_work;
                                        *spgiter = spgiter_work;
                                        *spgfcnt = spgfcnt_work;
                                        *tniter = tniter_work;
                                        *tnfcnt = tnfcnt_work;
                                        *tnstpcnt = tnstpcnt_work;
                                        *tnintcnt = tnintcnt_work;
                                        *tnexgcnt = tnexgcnt_work;
                                        *tnexbcnt = tnexbcnt_work;
                                        *tnintfe = tnintfe_work;
                                        *tnexgfe = tnexgfe_work;
                                        *tnexbfe = tnexbfe_work;
                                        *inform = spg_line_inform;
                                        return;
                                    }

                                    int grad_after_retry = 0;
                                    eval_gradient_full_cpp(
                                        n, x_work.data(), m, lambda, rho, gtype, g_work.data(), sterel, steabs, &grad_after_retry
                                    );
                                    gcnt_work += 1;

                                    if (grad_after_retry < 0) {
                                        *f = f_work;
                                        for (int i = 0; i < n_val; ++i) {
                                            x[i] = x_work[i];
                                            g[i] = g_work[i];
                                        }
                                        *iter = iter_work;
                                        *fcnt = fcnt_work;
                                        *gcnt = gcnt_work;
                                        *cgcnt = cgcnt_work;
                                        *spgiter = spgiter_work;
                                        *spgfcnt = spgfcnt_work;
                                        *tniter = tniter_work;
                                        *tnfcnt = tnfcnt_work;
                                        *tnstpcnt = tnstpcnt_work;
                                        *tnintcnt = tnintcnt_work;
                                        *tnexgcnt = tnexgcnt_work;
                                        *tnexbcnt = tnexbcnt_work;
                                        *tnintfe = tnintfe_work;
                                        *tnexgfe = tnexgfe_work;
                                        *tnexbfe = tnexbfe_work;
                                        *inform = grad_after_retry;
                                        return;
                                    }

                                    gpsupn_after = 0.0;
                                    gpeucn2_after = 0.0;
                                    for (int i = 0; i < n_val; ++i) {
                                        if (x_work[i] <= l_work[i] + std::max((*epsrel) * std::abs(l_work[i]), *epsabs)) {
                                            x_work[i] = l_work[i];
                                        } else if (x_work[i] >= u_work[i] - std::max((*epsrel) * std::abs(u_work[i]), *epsabs)) {
                                            x_work[i] = u_work[i];
                                        }
                                        const double xpg = x_work[i] - g_work[i];
                                        const double gpi = std::min(u_work[i], std::max(l_work[i], xpg)) - x_work[i];
                                        gpsupn_after = std::max(gpsupn_after, std::abs(gpi));
                                        gpeucn2_after += gpi * gpi;
                                    }

                                    int precision_after_retry_flag = 0;
                                    const bool precision_after_retry = packmolprecision_cpp(
                                        n, x_work.data(), m, lambda, rho, &precision_after_retry_flag
                                    );
                                    if (precision_after_retry_flag < 0) {
                                        *f = f_work;
                                        for (int i = 0; i < n_val; ++i) {
                                            x[i] = x_work[i];
                                            g[i] = g_work[i];
                                        }
                                        *iter = iter_work;
                                        *fcnt = fcnt_work;
                                        *gcnt = gcnt_work;
                                        *cgcnt = cgcnt_work;
                                        *spgiter = spgiter_work;
                                        *spgfcnt = spgfcnt_work;
                                        *tniter = tniter_work;
                                        *tnfcnt = tnfcnt_work;
                                        *tnstpcnt = tnstpcnt_work;
                                        *tnintcnt = tnintcnt_work;
                                        *tnexgcnt = tnexgcnt_work;
                                        *tnexbcnt = tnexbcnt_work;
                                        *tnintfe = tnintfe_work;
                                        *tnexgfe = tnexgfe_work;
                                        *tnexbfe = tnexbfe_work;
                                        *inform = precision_after_retry_flag;
                                        return;
                                    }
                                    post_inform = evaluate_post_step_inform_cpp(
                                        precision_after_retry,
                                        spg_line_inform,
                                        gpeucn2_after,
                                        *epsgpen,
                                        gpsupn_after,
                                        *epsgpsn,
                                        f_before_retry,
                                        f_work,
                                        *epsnfp,
                                        *maxitnfp,
                                        *infabs,
                                        lastgpns,
                                        *maxitngp,
                                        *fmin,
                                        iter_work,
                                        *maxit,
                                        fcnt_work,
                                        *maxfc,
                                        &progress_fprev,
                                        &progress_bestprog,
                                        &progress_itnfp
                                    );

                                    if (gencan_debug_enabled()) {
                                        std::fprintf(
                                            stderr,
                                            "[gencan-cpp-tn-retry] step_left=%d line_inform=%d post_inform=%d f=%.16e gpsupn=%.16e gpeucn2=%.16e\n",
                                            retry_budget,
                                            spg_line_inform,
                                            post_inform,
                                            f_work,
                                            gpsupn_after,
                                            gpeucn2_after
                                        );
                                    }
                                }

                                tn_post_debug_captured = true;
                                tn_post_line_inform = line_inform;
                                tn_post_post_inform = post_inform;
                                tn_post_nind = nind_work;
                                tn_post_f_before = f_try;
                                tn_post_f_after = f_work;
                                tn_post_gpsupn = gpsupn_after;
                                tn_post_gpeucn2 = gpeucn2_after;

                                if (post_inform >= 0) {
                                    *f = f_work;
                                    for (int i = 0; i < n_val; ++i) {
                                        x[i] = x_work[i];
                                        g[i] = g_work[i];
                                        s[i] = 0.0;
                                        y[i] = 0.0;
                                        d[i] = 0.0;
                                    }
                                    int nind_after = 0;
                                    for (int i = 0; i < n_val; ++i) {
                                        if (x_work[i] > l[i] && x_work[i] < u[i]) {
                                            ind[nind_after] = i + 1;
                                            nind_after += 1;
                                        }
                                    }
                                    *gpeucn2 = gpeucn2_after;
                                    *gpsupn = gpsupn_after;
                                    *iter = iter_work;
                                    *fcnt = fcnt_work;
                                    *gcnt = gcnt_work;
                                    *cgcnt = cgcnt_work;
                                    *spgiter = spgiter_work;
                                    *spgfcnt = spgfcnt_work;
                                    *tniter = tniter_work;
                                    *tnfcnt = tnfcnt_work;
                                    *tnstpcnt = tnstpcnt_work;
                                    *tnintcnt = tnintcnt_work;
                                    *tnexgcnt = tnexgcnt_work;
                                    *tnexbcnt = tnexbcnt_work;
                                    *tnintfe = tnintfe_work;
                                    *tnexgfe = tnexgfe_work;
                                    *tnexbfe = tnexbfe_work;
                                    if (*maxitngp > 0) {
                                        lastgpns[iter_work % (*maxitngp)] = gpeucn2_after;
                                    }
                                    *inform = post_inform;
                                    return;
                                }
                                fallback_reason = "tn_post_nonterminal";
                                const bool nested_tn_replay_seed =
                                    g_spg_post_cpp_replay_depth > 0 &&
                                    g_tn_post_cpp_replay_depth < tn_post_cpp_replay_max_depth();
                                const std::vector<double>& tn_x_seed = nested_tn_replay_seed ? x_try : x_work;
                                const std::vector<double>& tn_g_seed = nested_tn_replay_seed ? g_try : g_work;
                                seed_fallback_state(
                                    "tn_post_nonterminal",
                                    tn_x_seed,
                                    &tn_g_seed,
                                    nested_tn_replay_seed ? f_try : f_work
                                );
                                seed_fallback_cg_schedule(
                                    gencan_cg_schedule_seed_valid,
                                    gencan_cg_schedule_seed_valid ? gencan_cg_schedule_acgeps : acgeps,
                                    gencan_cg_schedule_seed_valid ? gencan_cg_schedule_bcgeps : bcgeps,
                                    gencan_cg_schedule_seed_valid ? gencan_cg_schedule_gpeucn20 : gpeucn20,
                                    gencan_cg_schedule_seed_valid ? gencan_cg_schedule_gpsupn0 : gpsupn0
                                );
                                seed_fallback_counters(
                                    nested_tn_replay_seed ? pre_tn_iter_seed : iter_work,
                                    nested_tn_replay_seed ? pre_tn_fcnt_seed : fcnt_work,
                                    nested_tn_replay_seed ? pre_tn_gcnt_seed : gcnt_work,
                                    nested_tn_replay_seed ? pre_tn_cgcnt_seed : cgcnt_work,
                                    nested_tn_replay_seed ? pre_tn_spgiter_seed : spgiter_work,
                                    nested_tn_replay_seed ? pre_tn_spgfcnt_seed : spgfcnt_work,
                                    nested_tn_replay_seed ? pre_tn_tniter_seed : tniter_work,
                                    nested_tn_replay_seed ? pre_tn_tnfcnt_seed : tnfcnt_work,
                                    nested_tn_replay_seed ? pre_tn_tnstpcnt_seed : tnstpcnt_work,
                                    nested_tn_replay_seed ? pre_tn_tnintcnt_seed : tnintcnt_work,
                                    nested_tn_replay_seed ? pre_tn_tnexgcnt_seed : tnexgcnt_work,
                                    nested_tn_replay_seed ? pre_tn_tnexbcnt_seed : tnexbcnt_work,
                                    nested_tn_replay_seed ? pre_tn_tnintfe_seed : tnintfe_work,
                                    nested_tn_replay_seed ? pre_tn_tnexgfe_seed : tnexgfe_work,
                                    nested_tn_replay_seed ? pre_tn_tnexbfe_seed : tnexbfe_work
                                );
                                seed_fallback_zero_steps();
                                if (!nested_tn_replay_seed) {
                                    for (int i = 0; i < n_val; ++i) {
                                        fallback_s_seed[i] = x_work[i] - x_try[i];
                                        fallback_y_seed[i] = g_work[i] - g_try[i];
                                    }
                                }
                                seed_fallback_progress(
                                    nested_tn_replay_seed ? pre_tn_progress_fprev : progress_fprev,
                                    nested_tn_replay_seed ? pre_tn_progress_bestprog : progress_bestprog,
                                    nested_tn_replay_seed ? pre_tn_progress_itnfp : progress_itnfp
                                );
                            }
                        }
                    } else {
                        fallback_reason = "tn_no_free_variables";
                        seed_fallback_state("tn_no_free_variables", x_try, &g_try, f_try);
                        seed_fallback_counters(
                            *iter, *fcnt, *gcnt, *cgcnt, *spgiter, *spgfcnt, *tniter, *tnfcnt,
                            *tnstpcnt, *tnintcnt, *tnexgcnt, *tnexbcnt, *tnintfe, *tnexgfe, *tnexbfe
                        );
                        seed_fallback_zero_steps();
                        seed_fallback_progress(progress_fprev, progress_bestprog, progress_itnfp);
                    }
                }
            }

            if (fallback_reason == nullptr) {
                if (gencan_debug_enabled()) {
                    std::fprintf(
                        stderr,
                        "[gencan-cpp-missing-tail-reason] mode=%d\n",
                        static_cast<int>(active_impl_mode())
                    );
                }
                *inform = 6;
                return;
            }

            const TailFallbackPolicy tail_policy = make_tail_fallback_policy(fallback_reason);
            const FallbackPostDebugInfo tn_post_debug{
                tn_post_debug_captured,
                tn_post_line_inform,
                tn_post_post_inform,
                tn_post_nind,
                tn_post_f_before,
                tn_post_f_after,
                tn_post_gpsupn,
                tn_post_gpeucn2
            };
            const FallbackPostDebugInfo spg_post_debug{
                spg_post_debug_captured,
                spg_post_line_inform,
                spg_post_post_inform,
                spg_post_nind,
                spg_post_f_before,
                spg_post_f_after,
                spg_post_gpsupn,
                spg_post_gpeucn2
            };
            emit_fallback_debug_cpp(tail_policy, fallback_reason, tn_post_debug, spg_post_debug);
            const bool force_blocked_seed_state = tail_policy.effective_block_tail;
            maybe_apply_fallback_seed_state_cpp(
                fallback_seed_state_enabled() || force_blocked_seed_state,
                fallback_x_seed_valid,
                n_val,
                fallback_x_seed,
                fallback_x_seed_reason,
                x
            );
            maybe_apply_fallback_gradient_seed_state_cpp(
                fallback_seed_state_enabled() || force_blocked_seed_state,
                fallback_g_seed_valid,
                n_val,
                fallback_g_seed,
                fallback_g_seed_reason,
                g
            );
            maybe_apply_fallback_f_seed_state_cpp(
                force_blocked_seed_state,
                fallback_f_seed_valid,
                fallback_f_seed,
                fallback_f_seed_reason,
                f
            );
            if (tail_policy.effective_block_tail && fallback_counter_seed_valid) {
                *iter = fallback_iter_seed;
                *fcnt = fallback_fcnt_seed;
                *gcnt = fallback_gcnt_seed;
                *cgcnt = fallback_cgcnt_seed;
                *spgiter = fallback_spgiter_seed;
                *spgfcnt = fallback_spgfcnt_seed;
                *tniter = fallback_tniter_seed;
                *tnfcnt = fallback_tnfcnt_seed;
                *tnstpcnt = fallback_tnstpcnt_seed;
                *tnintcnt = fallback_tnintcnt_seed;
                *tnexgcnt = fallback_tnexgcnt_seed;
                *tnexbcnt = fallback_tnexbcnt_seed;
                *tnintfe = fallback_tnintfe_seed;
                *tnexgfe = fallback_tnexgfe_seed;
                *tnexbfe = fallback_tnexbfe_seed;
            }
            if (tail_policy.effective_block_tail && fallback_step_seed_valid) {
                for (int i = 0; i < n_val; ++i) {
                    s[i] = (i < static_cast<int>(fallback_s_seed.size())) ? fallback_s_seed[i] : 0.0;
                    y[i] = (i < static_cast<int>(fallback_y_seed.size())) ? fallback_y_seed[i] : 0.0;
                }
            }
            const TnShadowSnapshot tn_shadow = capture_tn_shadow_snapshot_cpp(
                tn_post_shadow_enabled() && std::string(fallback_reason) == "tn_post_nonterminal",
                n_val,
                x,
                g,
                f,
                inform,
                iter,
                fcnt,
                gcnt,
                cgcnt,
                spgiter,
                tniter
            );
            const std::string reason(fallback_reason);
            if (reason == "spg_post_nonterminal") {
                continue_spg_post_nonterminal_cpp(
                    n_val,
                    n,
                    x,
                    l,
                    u,
                    m,
                    lambda,
                    rho,
                    epsgpen,
                    epsgpsn,
                    maxitnfp,
                    epsnfp,
                    maxitngp,
                    fmin,
                    maxit,
                    maxfc,
                    udelta0,
                    ucgmaxit,
                    cgscre,
                    cggpnf,
                    cgepsi,
                    cgepsf,
                    epsnqmp,
                    maxitnqmp,
                    nearlyq,
                    nint,
                    next,
                    mininterp,
                    maxextrap,
                    gtype,
                    htvtype,
                    trtype,
                    iprint,
                    ncomp,
                    f,
                    g,
                    gpeucn2,
                    gpsupn,
                    iter,
                    fcnt,
                    gcnt,
                    cgcnt,
                    spgiter,
                    spgfcnt,
                    tniter,
                    tnfcnt,
                    tnstpcnt,
                    tnintcnt,
                    tnexgcnt,
                    tnexbcnt,
                    tnintfe,
                    tnexgfe,
                    tnexbfe,
                    inform,
                    s,
                    y,
                    d,
                    ind,
                    lastgpns,
                    w,
                    eta,
                    delmin,
                    lspgma,
                    lspgmi,
                    theta,
                    gamma,
                    beta,
                    sigma1,
                    sigma2,
                    sterel,
                    steabs,
                    epsrel,
                    epsabs,
                    infrel,
                    infabs,
                    fallback_progress_fprev,
                    fallback_progress_bestprog,
                    fallback_progress_itnfp
                );
            } else if (is_tn_tail_reason_cpp(reason)) {
                continue_tn_tail_reason_cpp(
                    fallback_reason, n_val, tn_shadow,
                    fallback_progress_seed_valid, fallback_progress_fprev, fallback_progress_bestprog,
                    fallback_progress_itnfp, fallback_cg_schedule_seed_valid, fallback_cg_schedule_acgeps,
                    fallback_cg_schedule_bcgeps, fallback_cg_schedule_gpeucn20, fallback_cg_schedule_gpsupn0,
                    n, x, l, u, m, lambda, rho,
                    epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit, maxfc, udelta0,
                    ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq, nint,
                    next, mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g,
                    gpeucn2, gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt,
                    tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform, s,
                    y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta, gamma, beta, sigma1,
                    sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
                );
            } else {
                std::fprintf(
                    stderr,
                    "[gencan-cpp-missing-tail-reason] mode=%d\n",
                    static_cast<int>(active_impl_mode())
                );
            }
            return;
        }
    }
}
