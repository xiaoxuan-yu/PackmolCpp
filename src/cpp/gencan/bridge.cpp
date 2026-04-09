#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "gencan/internal.hpp"

thread_local int g_evalhd_fortran_call_count = 0;

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

bool ab_compare_enabled() {
    const char* env = std::getenv("PACKMOL_GENCAN_AB_COMPARE");
    if (env == nullptr) {
        return false;
    }
    return env[0] == '1' || env[0] == 't' || env[0] == 'T' || env[0] == 'y' || env[0] == 'Y';
}

bool is_tn_tail_reason_cpp(const std::string& reason) {
    return reason == "tn_post_nonterminal" || reason == "tn_no_free_variables";
}

enum class PendingContinuationPhase {
    kNone = 0,
    kSpgPost,
    kTnPost,
    kTnNoFreeVariables,
};

const char* pending_continuation_phase_reason_cpp(const PendingContinuationPhase phase) {
    switch (phase) {
    case PendingContinuationPhase::kSpgPost:
        return "spg_post_nonterminal";
    case PendingContinuationPhase::kTnPost:
        return "tn_post_nonterminal";
    case PendingContinuationPhase::kTnNoFreeVariables:
        return "tn_no_free_variables";
    case PendingContinuationPhase::kNone:
    default:
        return nullptr;
    }
}

constexpr int kSpgPostCppReplayMaxDepth = 1;
constexpr int kTnPostCppReplayMaxDepth = 1;

void apply_spg_post_cpp_replay_bookkeeping(
    int* iter,
    int* fcnt,
    int* gcnt
) {
    // Replay re-enters from the post-line-search state. The recursive bridge
    // already accounts for the continuation iteration itself, so the outer
    // handoff only needs to restore the pre-replay evaluation window.
    *fcnt += 3;
    *gcnt += 2;
}

double norm2_kernel(const int n, const double* x);

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
    const PendingContinuationPhase phase,
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
    if (phase != PendingContinuationPhase::kSpgPost) {
        return false;
    }
    if (active_impl_mode() != GencanImplMode::kCpp || !ab_compare_enabled()) {
        return false;
    }
    const int max_depth = kSpgPostCppReplayMaxDepth;
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
    const PendingContinuationPhase phase,
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
    const char* fallback_reason = pending_continuation_phase_reason_cpp(phase);
    const std::string reason = fallback_reason == nullptr ? "" : std::string(fallback_reason);
    if (phase == PendingContinuationPhase::kNone || !is_tn_tail_reason_cpp(reason)) {
        return false;
    }
    const bool nested_spg_replay = g_spg_post_cpp_replay_depth > 0;
    if (nested_spg_replay) {
        return false;
    }
    const bool cpp_mode = active_impl_mode() == GencanImplMode::kCpp && !ab_compare_enabled();
    const int max_depth = kTnPostCppReplayMaxDepth;
    const bool depth_ok = g_tn_post_cpp_replay_depth < max_depth;
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp-tn-post-cpp-replay-gate] reason=%s mode=%d replay_enabled=%d cpp_mode=%d depth=%d max_depth=%d\n",
            reason.c_str(),
            static_cast<int>(active_impl_mode()),
            cpp_mode ? 1 : 0,
            cpp_mode ? 1 : 0,
            g_tn_post_cpp_replay_depth,
            max_depth
        );
    }
    if (!cpp_mode || !depth_ok) {
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

extern "C" void packmol_gencan_reset_evalhd_call_count_c() {
    g_evalhd_fortran_call_count = 0;
}

extern "C" int packmol_gencan_evalhd_call_count_c() {
    return g_evalhd_fortran_call_count;
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

extern "C" void packmol_evalhd_fortran_c(
    const int* n
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
            bool gencan_cg_schedule_seed_valid = false;
            double gencan_cg_schedule_acgeps = 0.0;
            double gencan_cg_schedule_bcgeps = 0.0;
            double gencan_cg_schedule_gpeucn20 = 0.0;
            double gencan_cg_schedule_gpsupn0 = 0.0;
            auto apply_seeded_counters = [&](const int iter_seed, const int fcnt_seed, const int gcnt_seed,
                                             const int cgcnt_seed, const int spgiter_seed,
                                             const int spgfcnt_seed, const int tniter_seed,
                                             const int tnfcnt_seed, const int tnstpcnt_seed,
                                             const int tnintcnt_seed, const int tnexgcnt_seed,
                                             const int tnexbcnt_seed, const int tnintfe_seed,
                                             const int tnexgfe_seed, const int tnexbfe_seed) {
                *iter = iter_seed;
                *fcnt = fcnt_seed;
                *gcnt = gcnt_seed;
                *cgcnt = cgcnt_seed;
                *spgiter = spgiter_seed;
                *spgfcnt = spgfcnt_seed;
                *tniter = tniter_seed;
                *tnfcnt = tnfcnt_seed;
                *tnstpcnt = tnstpcnt_seed;
                *tnintcnt = tnintcnt_seed;
                *tnexgcnt = tnexgcnt_seed;
                *tnexbcnt = tnexbcnt_seed;
                *tnintfe = tnintfe_seed;
                *tnexgfe = tnexgfe_seed;
                *tnexbfe = tnexbfe_seed;
            };
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

                    if (ab_compare_enabled() && try_spg_post_cpp_replay_tail_cpp(
                                   PendingContinuationPhase::kSpgPost, n, x, l, u, m, lambda, rho, epsgpen,
                                   epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit, maxfc,
                                   udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp,
                                   maxitnqmp, nearlyq, nint, next, mininterp, maxextrap, gtype,
                                   htvtype, trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter,
                                   fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt,
                                   tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe,
                                   inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma,
                                   lspgmi, theta, gamma, beta, sigma1, sigma2, sterel, steabs,
                                   epsrel, epsabs, infrel, infabs
                               )) {
                        return;
                    }
                    *inform = finalize_post_step_inform_cpp(
                        post_inform,
                        precision_after_spg,
                        gpeucn2_after,
                        *epsgpen,
                        gpsupn_after,
                        *epsgpsn,
                        iter_work,
                        *maxit,
                        fcnt_work,
                        *maxfc
                    );
                    return;
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
                                const bool nested_tn_replay_seed =
                                    g_spg_post_cpp_replay_depth > 0 &&
                                    g_tn_post_cpp_replay_depth < kTnPostCppReplayMaxDepth;
                                const std::vector<double>& tn_x_seed = nested_tn_replay_seed ? x_try : x_work;
                                const std::vector<double>& tn_g_seed = nested_tn_replay_seed ? g_try : g_work;
                                *f = nested_tn_replay_seed ? f_try : f_work;
                                for (int i = 0; i < n_val; ++i) {
                                    x[i] = tn_x_seed[i];
                                    g[i] = tn_g_seed[i];
                                    s[i] = 0.0;
                                    y[i] = 0.0;
                                }
                                apply_seeded_counters(
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
                                if (!nested_tn_replay_seed) {
                                    for (int i = 0; i < n_val; ++i) {
                                        s[i] = x_work[i] - x_try[i];
                                        y[i] = g_work[i] - g_try[i];
                                    }
                                }
                                if (try_tn_post_cpp_replay_tail_cpp(
                                        PendingContinuationPhase::kTnPost, n, x, l, u, m, lambda, rho,
                                        epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit, maxfc,
                                        udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp,
                                        maxitnqmp, nearlyq, nint, next, mininterp, maxextrap, gtype,
                                        htvtype, trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter,
                                        fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt,
                                        tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform,
                                        s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta,
                                        gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs,
                                        infrel, infabs
                                    )) {
                                } else {
                                    *inform = finalize_post_step_inform_cpp(
                                        post_inform,
                                        precision_after,
                                        gpeucn2_after,
                                        *epsgpen,
                                        gpsupn_after,
                                        *epsgpsn,
                                        iter_work,
                                        *maxit,
                                        fcnt_work,
                                        *maxfc
                                    );
                                }
                                return;
                            }
                        }
                    } else {
                        *f = f_try;
                        for (int i = 0; i < n_val; ++i) {
                            x[i] = x_try[i];
                            g[i] = g_try[i];
                            s[i] = 0.0;
                            y[i] = 0.0;
                        }
                        for (int i = 0; i < nind_try; ++i) {
                            ind[i] = ind_try[i];
                        }
                        *gpeucn2 = gpeucn2_try;
                        *gpsupn = gpsupn_try;
                        if (*maxitngp > 0) {
                            lastgpns[*iter % (*maxitngp)] = gpeucn2_try;
                        }
                        apply_seeded_counters(
                            *iter, *fcnt, *gcnt, *cgcnt, *spgiter, *spgfcnt, *tniter, *tnfcnt,
                            *tnstpcnt, *tnintcnt, *tnexgcnt, *tnexbcnt, *tnintfe, *tnexgfe, *tnexbfe
                        );
                        *inform = finalize_post_step_inform_cpp(
                            -1,
                            precision_solution,
                            gpeucn2_try,
                            *epsgpen,
                            gpsupn_try,
                            *epsgpsn,
                            *iter,
                            *maxit,
                            *fcnt,
                            *maxfc
                        );
                        return;
                    }
                }
            }

        }
    }
}
