#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "gencan/api.hpp"

thread_local int g_evalhd_fortran_call_count = 0;

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

void packmol_gencan_tnls_bridge_impl(
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

void packmol_gencan_cg_bridge_impl(
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

namespace {

constexpr int kContinueMainLoop = 9;

struct SolverWorkCounters {
    int iter;
    int fcnt;
    int gcnt;
    int cgcnt;
    int spgiter;
    int spgfcnt;
    int tniter;
    int tnfcnt;
    int tnstpcnt;
    int tnintcnt;
    int tnexgcnt;
    int tnexbcnt;
    int tnintfe;
    int tnexgfe;
    int tnexbfe;
};

struct ProjectedGradientSummary {
    double gpsupn;
    double gpeucn2;
    double gieucn2;
    int nind;
};

struct TnReplayBuffers {
    std::vector<double> xplus;
    std::vector<double> xtmp;
    std::vector<double> xbext;
};

struct SolverCommitView {
    int n;
    GencanCounters counters;
    double* f;
    double* x;
    double* g;
    double* gpeucn2;
    double* gpsupn;
    int* inform;
    double* s;
    double* y;
    double* d;
    int* ind;
    double* lastgpns;
    double* w;
};

struct CgScheduleSeed {
    bool valid;
    double acgeps;
    double bcgeps;
    double gpeucn20;
    double gpsupn0;
};

SolverWorkCounters make_solver_work_counters(
    const int iter_seed,
    const int fcnt_seed,
    const int gcnt_seed,
    const int cgcnt_seed,
    const int spgiter_seed,
    const int spgfcnt_seed,
    const int tniter_seed,
    const int tnfcnt_seed,
    const int tnstpcnt_seed,
    const int tnintcnt_seed,
    const int tnexgcnt_seed,
    const int tnexbcnt_seed,
    const int tnintfe_seed,
    const int tnexgfe_seed,
    const int tnexbfe_seed
) {
    return SolverWorkCounters{
        iter_seed, fcnt_seed, gcnt_seed, cgcnt_seed, spgiter_seed, spgfcnt_seed, tniter_seed,
        tnfcnt_seed, tnstpcnt_seed, tnintcnt_seed, tnexgcnt_seed, tnexbcnt_seed, tnintfe_seed,
        tnexgfe_seed, tnexbfe_seed
    };
}

SolverWorkCounters load_solver_work_counters(const GencanCounters& counters) {
    return make_solver_work_counters(
        *counters.iter,
        *counters.fcnt,
        *counters.gcnt,
        *counters.cgcnt,
        *counters.spgiter,
        *counters.spgfcnt,
        *counters.tniter,
        *counters.tnfcnt,
        *counters.tnstpcnt,
        *counters.tnintcnt,
        *counters.tnexgcnt,
        *counters.tnexbcnt,
        *counters.tnintfe,
        *counters.tnexgfe,
        *counters.tnexbfe
    );
}

void clamp_to_box_inplace(
    const int n,
    double* x,
    const double* l,
    const double* u,
    const double epsrel,
    const double epsabs
) {
    for (int i = 0; i < n; ++i) {
        if (x[i] <= l[i] + std::max(epsrel * std::abs(l[i]), epsabs)) {
            x[i] = l[i];
        } else if (x[i] >= u[i] - std::max(epsrel * std::abs(u[i]), epsabs)) {
            x[i] = u[i];
        }
    }
}

ProjectedGradientSummary analyze_projected_gradient(
    const int n,
    const double* x,
    const double* g,
    const double* l,
    const double* u,
    int* ind_out
) {
    ProjectedGradientSummary summary{0.0, 0.0, 0.0, 0};
    for (int i = 0; i < n; ++i) {
        const double xpg = x[i] - g[i];
        const double gpi = std::min(u[i], std::max(l[i], xpg)) - x[i];
        summary.gpsupn = std::max(summary.gpsupn, std::abs(gpi));
        summary.gpeucn2 += gpi * gpi;
        if (x[i] > l[i] && x[i] < u[i]) {
            summary.gieucn2 += gpi * gpi;
            if (ind_out != nullptr) {
                ind_out[summary.nind] = i + 1;
            }
            summary.nind += 1;
        }
    }
    return summary;
}

void trace_solver_line(
    const char* tag,
    const double f_value,
    const PackmolPrecisionState& precision_state,
    const ProjectedGradientSummary& summary,
    const double xnorm,
    const double ometa2
) {
    if (!gencan_debug_enabled()) {
        return;
    }
    std::fprintf(
        stderr,
        "[gencan-cpp] %s f=%.17e precision=%.17e fdist=%.17e frest=%.17e gpsupn=%.17e gpeucn2=%.17e gieucn2=%.17e xnorm=%.17e spg_test=%d\n",
        tag,
        f_value,
        precision_state.precision_value,
        precision_state.fdist_value,
        precision_state.frest_value,
        summary.gpsupn,
        summary.gpeucn2,
        summary.gieucn2,
        xnorm,
        (summary.gpeucn2 > 0.0 && summary.gieucn2 <= ometa2 * summary.gpeucn2) ? 1 : 0
    );
}

void apply_solver_work_counters(
    const GencanCounters& target,
    const SolverWorkCounters& work_counters
) {
    *target.iter = work_counters.iter;
    *target.fcnt = work_counters.fcnt;
    *target.gcnt = work_counters.gcnt;
    *target.cgcnt = work_counters.cgcnt;
    *target.spgiter = work_counters.spgiter;
    *target.spgfcnt = work_counters.spgfcnt;
    *target.tniter = work_counters.tniter;
    *target.tnfcnt = work_counters.tnfcnt;
    *target.tnstpcnt = work_counters.tnstpcnt;
    *target.tnintcnt = work_counters.tnintcnt;
    *target.tnexgcnt = work_counters.tnexgcnt;
    *target.tnexbcnt = work_counters.tnexbcnt;
    *target.tnintfe = work_counters.tnintfe;
    *target.tnexgfe = work_counters.tnexgfe;
    *target.tnexbfe = work_counters.tnexbfe;
}

void commit_xy_state(
    const SolverCommitView& view,
    const double f_value,
    const std::vector<double>& x_src,
    const std::vector<double>& g_src,
    const SolverWorkCounters& work_counters,
    const int inform_value
) {
    *view.f = f_value;
    for (int i = 0; i < view.n; ++i) {
        view.x[i] = x_src[i];
        view.g[i] = g_src[i];
    }
    apply_solver_work_counters(view.counters, work_counters);
    *view.inform = inform_value;
}

void commit_full_state(
    const SolverCommitView& view,
    const double f_value,
    const std::vector<double>& x_src,
    const std::vector<double>& g_src,
    const std::vector<double>& s_src,
    const std::vector<double>& y_src,
    const std::vector<double>& d_src,
    const std::vector<int>& ind_src,
    const int nind_src,
    const double gpeucn2_value,
    const double gpsupn_value,
    const SolverWorkCounters& work_counters,
    const int inform_value
) {
    *view.f = f_value;
    for (int i = 0; i < view.n; ++i) {
        view.x[i] = x_src[i];
        view.g[i] = g_src[i];
        view.s[i] = s_src[i];
        view.y[i] = y_src[i];
        view.d[i] = d_src[i];
    }
    for (int i = 0; i < nind_src; ++i) {
        view.ind[i] = ind_src[i];
    }
    *view.gpeucn2 = gpeucn2_value;
    *view.gpsupn = gpsupn_value;
    apply_solver_work_counters(view.counters, work_counters);
    *view.inform = inform_value;
}

double compute_lamspg(
    const double xnorm,
    const double gpeucn2_value,
    const GencanOptions& options
) {
    double lamspg = std::max(1.0, xnorm) / std::sqrt(gpeucn2_value);
    lamspg = std::min(*options.lspgma, std::max(*options.lspgmi, lamspg));
    return lamspg;
}

int run_spg_post_phase(
    const GencanProblemView& problem,
    const GencanOptions& options,
    const SolverCommitView& commit_view,
    const int iter_value,
    const double xnorm_before,
    const double sts_before,
    const double sty_before,
    const double f_try,
    const double gpeucn2_try,
    const std::vector<double>& x_try,
    const std::vector<double>& g_try
) {
    std::vector<double> x_work = x_try;
    std::vector<double> g_work = g_try;
    std::vector<double> x_prev = x_try;
    std::vector<double> g_prev = g_try;
    std::vector<double> s_work(commit_view.n, 0.0);
    std::vector<double> y_work(commit_view.n, 0.0);
    std::vector<double> xtrial_work(commit_view.n, 0.0);
    std::vector<double> d_work(commit_view.n, 0.0);
    std::vector<int> ind_work(commit_view.n, 0);

    SolverWorkCounters spg_counters = load_solver_work_counters(commit_view.counters);
    spg_counters.spgiter += 1;
    auto current_counters = [&]() { return spg_counters; };
    double f_work = f_try;

    double lamspg = 0.0;
    if (iter_value == 1 || sty_before <= 0.0) {
        lamspg = compute_lamspg(xnorm_before, gpeucn2_try, options);
    } else {
        lamspg = std::min(*options.lspgma, std::max(*options.lspgmi, sts_before / sty_before));
    }

    int ls_inform = 0;
    const int fcnt_prev = spg_counters.fcnt;
    run_spgls_context_cpp(GencanSpglsContext{
        problem.n, x_work.data(), problem.m, problem.lambda, problem.rho, &f_work, g_work.data(),
        problem.l, problem.u, &lamspg, options.nint, options.mininterp, options.fmin, options.maxfc,
        problem.iprint, &spg_counters.fcnt, &ls_inform, xtrial_work.data(), d_work.data(),
        options.gamma, options.sigma1, options.sigma2, options.sterel, options.steabs,
        options.epsrel, options.epsabs, options.infrel, options.infabs
    });
    spg_counters.spgfcnt += (spg_counters.fcnt - fcnt_prev);

    if (ls_inform < 0) {
        commit_xy_state(commit_view, f_work, x_work, g_work, current_counters(), ls_inform);
        return ls_inform;
    }

    int grad_after_spg = 0;
    eval_gradient_full_cpp(
        problem.n, x_work.data(), problem.m, problem.lambda, problem.rho, problem.gtype,
        g_work.data(), options.sterel, options.steabs, &grad_after_spg
    );
    spg_counters.gcnt += 1;

    if (grad_after_spg < 0) {
        commit_xy_state(commit_view, f_work, x_work, g_work, current_counters(), grad_after_spg);
        return grad_after_spg;
    }

    clamp_to_box_inplace(commit_view.n, x_work.data(), problem.l, problem.u, *options.epsrel, *options.epsabs);

    for (int i = 0; i < commit_view.n; ++i) {
        s_work[i] = x_work[i] - x_prev[i];
        y_work[i] = g_work[i] - g_prev[i];
    }
    const ProjectedGradientSummary spg_summary = analyze_projected_gradient(
        commit_view.n, x_work.data(), g_work.data(), problem.l, problem.u, ind_work.data()
    );
    const double gpsupn_after = spg_summary.gpsupn;
    const double gpeucn2_after = spg_summary.gpeucn2;
    const int nind_after = spg_summary.nind;

    if (ls_inform == 6) {
        commit_full_state(
            commit_view, f_work, x_work, g_work, s_work, y_work, d_work, ind_work, nind_after,
            gpeucn2_after, gpsupn_after, current_counters(), 6
        );
        return 6;
    }

    commit_full_state(
        commit_view, f_work, x_work, g_work, s_work, y_work, d_work, ind_work, nind_after,
        gpeucn2_after, gpsupn_after, current_counters(), -1
    );
    return kContinueMainLoop;
}

int run_tn_post_phase(
    const GencanProblemView& problem,
    const GencanOptions& options,
    const SolverCommitView& commit_view,
    const CgScheduleSeed& cg_schedule_seed,
    const int iter_value,
    const double xnorm_before,
    const double sts_before,
    const double sty_before,
    const double f_try,
    const double gpsupn_try,
    const double gpeucn2_try,
    const int nind_try,
    const std::vector<double>& x_try,
    const std::vector<double>& g_try,
    const std::vector<int>& ind_try
) {
    std::vector<double> x_work = x_try;
    std::vector<double> g_work = g_try;
    std::vector<double> l_work(commit_view.n);
    std::vector<double> u_work(commit_view.n);
    for (int i = 0; i < commit_view.n; ++i) {
        l_work[i] = problem.l[i];
        u_work[i] = problem.u[i];
    }

    const int nind_work = nind_try;
    if (nind_work <= 0) {
        *commit_view.f = f_try;
        for (int i = 0; i < commit_view.n; ++i) {
            commit_view.x[i] = x_try[i];
            commit_view.g[i] = g_try[i];
            commit_view.s[i] = 0.0;
            commit_view.y[i] = 0.0;
        }
        for (int i = 0; i < nind_try; ++i) {
            commit_view.ind[i] = ind_try[i];
        }
        *commit_view.gpeucn2 = gpeucn2_try;
        *commit_view.gpsupn = gpsupn_try;
        apply_solver_work_counters(commit_view.counters, load_solver_work_counters(commit_view.counters));
        *commit_view.inform = -1;
        return kContinueMainLoop;
    }

    shrink_inplace(nind_work, ind_try.data(), x_work.data());
    shrink_inplace(nind_work, ind_try.data(), g_work.data());
    shrink_inplace(nind_work, ind_try.data(), l_work.data());
    shrink_inplace(nind_work, ind_try.data(), u_work.data());

    SolverWorkCounters tn_work_counters = load_solver_work_counters(commit_view.counters);
    tn_work_counters.tniter += 1;
    int& fcnt_work = tn_work_counters.fcnt;
    int& gcnt_work = tn_work_counters.gcnt;
    int& cgcnt_work = tn_work_counters.cgcnt;
    int& spgiter_work = tn_work_counters.spgiter;
    int& spgfcnt_work = tn_work_counters.spgfcnt;
    int& tnfcnt_work = tn_work_counters.tnfcnt;
    int& tnstpcnt_work = tn_work_counters.tnstpcnt;
    int& tnintcnt_work = tn_work_counters.tnintcnt;
    int& tnexgcnt_work = tn_work_counters.tnexgcnt;
    int& tnexbcnt_work = tn_work_counters.tnexbcnt;
    int& tnintfe_work = tn_work_counters.tnintfe;
    int& tnexgfe_work = tn_work_counters.tnexgfe;
    int& tnexbfe_work = tn_work_counters.tnexbfe;
    double delta = 0.0;
    if (iter_value == 1) {
        delta = (*options.udelta0 <= 0.0)
            ? std::max(*options.delmin, 0.1 * std::max(1.0, xnorm_before))
            : *options.udelta0;
    } else {
        delta = std::max(*options.delmin, 10.0 * std::sqrt(sts_before));
    }

    const double epsgpen2 = (*options.epsgpen) * (*options.epsgpen);
    double acgeps = 0.0;
    double bcgeps = 0.0;
    if (cg_schedule_seed.valid) {
        acgeps = cg_schedule_seed.acgeps;
        bcgeps = cg_schedule_seed.bcgeps;
    } else {
        gp_ieee_signal1_cpp(
            gpsupn_try, &acgeps, &bcgeps, *options.cgepsf, *options.cgepsi, *options.cggpnf
        );
    }
    const double gpeucn20 = cg_schedule_seed.valid ? cg_schedule_seed.gpeucn20 : gpeucn2_try;
    const double gpsupn0 = cg_schedule_seed.valid ? cg_schedule_seed.gpsupn0 : gpsupn_try;
    double kappa = 0.0;
    double cgeps = *options.cgepsf;
    int cgmaxit = 0;
    gp_ieee_signal2_cpp(
        &cgmaxit, nind_work, *options.nearlyq, *options.ucgmaxit, *options.cgscre, &kappa, gpeucn2_try,
        gpeucn20, epsgpen2, *options.epsgpsn, &cgeps, acgeps, bcgeps, *options.cgepsf,
        *options.cgepsi, gpsupn_try, gpsupn0
    );
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp] tn-setup iter=%d delta=%.17e cgeps=%.17e cgmaxit=%d kappa=%.17e gpsupn=%.17e gpeucn2=%.17e xnorm=%.17e\n",
            iter_value,
            delta,
            cgeps,
            cgmaxit,
            kappa,
            gpsupn_try,
            gpeucn2_try,
            xnorm_before
        );
    }

    std::vector<double> cg_s(commit_view.n, 0.0);
    std::vector<double> cg_w(commit_view.n, 0.0);
    std::vector<double> cg_y(commit_view.n, 0.0);
    std::vector<double> cg_r(commit_view.n, 0.0);
    std::vector<double> cg_d(commit_view.n, 0.0);
    std::vector<double> cg_sprev(commit_view.n, 0.0);

    int cg_iter = 0;
    int rbdtype = 0;
    int rbdind = 0;
    int cg_inform = 0;
    run_cg_context_cpp(GencanCgContext{
        &nind_work, ind_try.data(), problem.n, x_work.data(), problem.m, problem.lambda, problem.rho,
        g_work.data(), &delta, l_work.data(), u_work.data(), &cgeps, options.epsnqmp,
        options.maxitnqmp, &cgmaxit, options.nearlyq, problem.gtype, problem.htvtype, problem.trtype,
        problem.iprint, problem.ncomp, cg_s.data(), &cg_iter, &rbdtype, &rbdind, &cg_inform,
        cg_w.data(), cg_y.data(), cg_r.data(), cg_d.data(), cg_sprev.data(), options.theta,
        options.sterel, options.steabs, options.epsrel, options.epsabs, options.infrel, options.infabs
    });
    cgcnt_work += cg_iter;
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp] tn-cg iter=%d inform=%d cg_iter=%d rbdtype=%d rbdind=%d\n",
            iter_value,
            cg_inform,
            cg_iter,
            rbdtype,
            rbdind
        );
    }

    if (cg_inform < 0) {
        commit_xy_state(commit_view, f_try, x_work, g_work, tn_work_counters, cg_inform);
        return cg_inform;
    }

    double amax = *options.infabs;
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
    TnReplayBuffers tn_buffers{
        std::vector<double>(commit_view.n, 0.0),
        std::vector<double>(commit_view.n, 0.0),
        std::vector<double>(commit_view.n, 0.0)
    };
    run_tnls_context_cpp(GencanTnlsContext{
        &nind_work, ind_try.data(), problem.n, x_work.data(), problem.m, problem.lambda, problem.rho,
        l_work.data(), u_work.data(), &f_work, g_work.data(), cg_s.data(), &amax, &rbdtype, &rbdind,
        options.nint, options.next, options.mininterp, options.maxextrap, options.fmin, options.maxfc,
        problem.gtype, problem.iprint, &fcnt_work, &gcnt_work, &tnintcnt_work, &tnexgcnt_work,
        &tnexbcnt_work, &tnls_inform, tn_buffers.xplus.data(), tn_buffers.xtmp.data(),
        tn_buffers.xbext.data(), options.gamma, options.beta, options.sigma1, options.sigma2,
        options.sterel, options.steabs, options.epsrel, options.epsabs, options.infrel, options.infabs
    });

    if (tnls_inform < 0) {
        commit_xy_state(commit_view, f_work, x_work, g_work, tn_work_counters, tnls_inform);
        return tnls_inform;
    }
    if (gencan_debug_enabled()) {
        std::fprintf(
            stderr,
            "[gencan-cpp] tn-line iter=%d inform=%d f=%.17e fcnt=%d gcnt=%d\n",
            iter_value,
            tnls_inform,
            f_work,
            fcnt_work,
            gcnt_work
        );
    }

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

    clamp_to_box_inplace(commit_view.n, x_work.data(), l_work.data(), u_work.data(), *options.epsrel, *options.epsabs);

    int line_inform = tnls_inform;
    if (tnls_inform == 6) {
        spgiter_work += 1;
        const double lamspg = (iter_value == 1 || sty_before <= 0.0)
            ? compute_lamspg(xnorm_before, gpeucn2_try, options)
            : std::min(*options.lspgma, std::max(*options.lspgmi, sts_before / sty_before));

        const int fcnt_prev_spg = fcnt_work;
        double lamspg_mut = lamspg;
        run_spgls_context_cpp(GencanSpglsContext{
            problem.n, x_work.data(), problem.m, problem.lambda, problem.rho, &f_work, g_work.data(),
            l_work.data(), u_work.data(), &lamspg_mut, options.nint, options.mininterp, options.fmin,
            options.maxfc, problem.iprint, &fcnt_work, &line_inform, tn_buffers.xplus.data(),
            tn_buffers.xtmp.data(), options.gamma, options.sigma1, options.sigma2, options.sterel,
            options.steabs, options.epsrel, options.epsabs, options.infrel, options.infabs
        });
        spgfcnt_work += (fcnt_work - fcnt_prev_spg);

        if (line_inform < 0) {
            commit_xy_state(commit_view, f_work, x_work, g_work, tn_work_counters, line_inform);
            return line_inform;
        }

        int grad_after_spg = 0;
        eval_gradient_full_cpp(
            problem.n, x_work.data(), problem.m, problem.lambda, problem.rho, problem.gtype,
            g_work.data(), options.sterel, options.steabs, &grad_after_spg
        );
        gcnt_work += 1;

        if (grad_after_spg < 0) {
            commit_xy_state(commit_view, f_work, x_work, g_work, tn_work_counters, grad_after_spg);
            return grad_after_spg;
        }
    }

    const ProjectedGradientSummary final_summary = analyze_projected_gradient(
        commit_view.n, x_work.data(), g_work.data(), problem.l, problem.u, commit_view.ind
    );
    *commit_view.f = f_work;
    for (int i = 0; i < commit_view.n; ++i) {
        commit_view.x[i] = x_work[i];
        commit_view.g[i] = g_work[i];
        commit_view.s[i] = x_work[i] - x_try[i];
        commit_view.y[i] = g_work[i] - g_try[i];
    }
    *commit_view.gpeucn2 = final_summary.gpeucn2;
    *commit_view.gpsupn = final_summary.gpsupn;
    apply_solver_work_counters(commit_view.counters, tn_work_counters);
    if (line_inform == 6) {
        *commit_view.inform = 6;
        return 6;
    }
    *commit_view.inform = -1;
    return kContinueMainLoop;
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

}  // namespace

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
) {
    const GencanProblemView problem{
        n, x, l, u, m, lambda, rho, gtype, htvtype, trtype, iprint, ncomp
    };
    const GencanOptions options{
        epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit, maxfc, udelta0, ucgmaxit,
        cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq, nint, next, mininterp,
        maxextrap, eta, delmin, lspgma, lspgmi, theta, gamma, beta, sigma1, sigma2, sterel,
        steabs, epsrel, epsabs, infrel, infabs
    };
    GencanCounters counters{
        iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt,
        tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe
    };
    GencanMutableState state{f, g, gpeucn2, gpsupn, inform, ind, lastgpns};
    GencanWorkspace workspace{s, y, d, w};

    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_gencan_fortran_c(
                problem.n, problem.x, problem.l, problem.u, problem.m, problem.lambda, problem.rho,
                options.epsgpen, options.epsgpsn, options.maxitnfp, options.epsnfp,
                options.maxitngp, options.fmin, options.maxit, options.maxfc, options.udelta0,
                options.ucgmaxit, options.cgscre, options.cggpnf, options.cgepsi, options.cgepsf,
                options.epsnqmp, options.maxitnqmp, options.nearlyq, options.nint, options.next,
                options.mininterp, options.maxextrap, problem.gtype, problem.htvtype, problem.trtype,
                problem.iprint, problem.ncomp, state.f, state.g, state.gpeucn2, state.gpsupn,
                counters.iter, counters.fcnt, counters.gcnt, counters.cgcnt, counters.spgiter,
                counters.spgfcnt, counters.tniter, counters.tnfcnt, counters.tnstpcnt,
                counters.tnintcnt, counters.tnexgcnt, counters.tnexbcnt, counters.tnintfe,
                counters.tnexgfe, counters.tnexbfe, state.inform, workspace.s, workspace.y,
                workspace.d, state.ind, state.lastgpns, workspace.w, options.eta, options.delmin,
                options.lspgma, options.lspgmi, options.theta, options.gamma, options.beta,
                options.sigma1, options.sigma2, options.sterel, options.steabs, options.epsrel,
                options.epsabs, options.infrel, options.infabs
            );
            return;
        case GencanImplMode::kCpp: {
            n = problem.n;
            x = problem.x;
            l = problem.l;
            u = problem.u;
            m = problem.m;
            lambda = problem.lambda;
            rho = problem.rho;
            gtype = problem.gtype;
            htvtype = problem.htvtype;
            trtype = problem.trtype;
            iprint = problem.iprint;
            ncomp = problem.ncomp;

            epsgpen = options.epsgpen;
            epsgpsn = options.epsgpsn;
            maxitnfp = options.maxitnfp;
            epsnfp = options.epsnfp;
            maxitngp = options.maxitngp;
            fmin = options.fmin;
            maxit = options.maxit;
            maxfc = options.maxfc;
            udelta0 = options.udelta0;
            ucgmaxit = options.ucgmaxit;
            cgscre = options.cgscre;
            cggpnf = options.cggpnf;
            cgepsi = options.cgepsi;
            cgepsf = options.cgepsf;
            epsnqmp = options.epsnqmp;
            maxitnqmp = options.maxitnqmp;
            nearlyq = options.nearlyq;
            nint = options.nint;
            next = options.next;
            mininterp = options.mininterp;
            maxextrap = options.maxextrap;
            eta = options.eta;
            delmin = options.delmin;
            lspgma = options.lspgma;
            lspgmi = options.lspgmi;
            theta = options.theta;
            gamma = options.gamma;
            beta = options.beta;
            sigma1 = options.sigma1;
            sigma2 = options.sigma2;
            sterel = options.sterel;
            steabs = options.steabs;
            epsrel = options.epsrel;
            epsabs = options.epsabs;
            infrel = options.infrel;
            infabs = options.infabs;

            f = state.f;
            g = state.g;
            gpeucn2 = state.gpeucn2;
            gpsupn = state.gpsupn;
            inform = state.inform;
            ind = state.ind;
            lastgpns = state.lastgpns;

            s = workspace.s;
            y = workspace.y;
            d = workspace.d;
            w = workspace.w;

            iter = counters.iter;
            fcnt = counters.fcnt;
            gcnt = counters.gcnt;
            cgcnt = counters.cgcnt;
            spgiter = counters.spgiter;
            spgfcnt = counters.spgfcnt;
            tniter = counters.tniter;
            tnfcnt = counters.tnfcnt;
            tnstpcnt = counters.tnstpcnt;
            tnintcnt = counters.tnintcnt;
            tnexgcnt = counters.tnexgcnt;
            tnexbcnt = counters.tnexbcnt;
            tnintfe = counters.tnintfe;
            tnexgfe = counters.tnexgfe;
            tnexbfe = counters.tnexbfe;

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
            CgScheduleSeed cg_schedule_seed{false, 0.0, 0.0, 0.0, 0.0};
            const SolverCommitView commit_view{
                n_val, counters, f, x, g, gpeucn2, gpsupn, inform, s, y, d, ind, lastgpns, w
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
            eval_objective_full_cpp(n, x_try.data(), m, lambda, rho, &f_try, &eval_flag);

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
            PackmolPrecisionState precision_state = read_packmol_precision_state_cpp();
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
                trace_solver_line(
                    "entry-precision-solution",
                    f_try,
                    precision_state,
                    ProjectedGradientSummary{0.0, 0.0, 0.0, 0},
                    std::sqrt(norm2_kernel(n_val, x_try.data())),
                    0.0
                );
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
                const ProjectedGradientSummary try_summary =
                    analyze_projected_gradient(n_val, x_try.data(), g_try.data(), l, u, ind_try.data());
                double gpsupn_try = try_summary.gpsupn;
                double gpeucn2_try = try_summary.gpeucn2;
                double gieucn2_try = try_summary.gieucn2;
                int nind_current = try_summary.nind;

                *f = f_try;
                for (int i = 0; i < n_val; ++i) {
                    x[i] = x_try[i];
                    g[i] = g_try[i];
                }
                for (int i = 0; i < nind_current; ++i) {
                    ind[i] = ind_try[i];
                }
                *gpeucn2 = gpeucn2_try;
                *gpsupn = gpsupn_try;

                const double epsgpen2 = (*epsgpen) * (*epsgpen);
                gp_ieee_signal1_cpp(
                    gpsupn_try,
                    &cg_schedule_seed.acgeps,
                    &cg_schedule_seed.bcgeps,
                    *cgepsf,
                    *cgepsi,
                    *cggpnf
                );
                cg_schedule_seed.gpeucn20 = gpeucn2_try;
                cg_schedule_seed.gpsupn0 = gpsupn_try;
                cg_schedule_seed.valid = true;

                if (*maxitngp > 0) {
                    for (int i = 0; i < *maxitngp; ++i) {
                        lastgpns[i] = *infabs;
                    }
                }

                double fprev = *infabs;
                double bestprog = 0.0;
                int itnfp = 0;
                double xnorm = std::sqrt(norm2_kernel(n_val, x));
                double sts = 0.0;
                double sty = 0.0;
                const double ometa2 = (1.0 - *eta) * (1.0 - *eta);
                trace_solver_line("entry", f_try, precision_state, try_summary, xnorm, ometa2);

                while (true) {
                    int precision_loop_flag = 0;
                    const bool precision_loop = packmolprecision_cpp(
                        n, x, m, lambda, rho, &precision_loop_flag
                    );
                    precision_state = read_packmol_precision_state_cpp();
                    if (precision_loop_flag < 0) {
                        *inform = precision_loop_flag;
                        return;
                    }
                    trace_solver_line("loop", *f, precision_state, ProjectedGradientSummary{*gpsupn, *gpeucn2, gieucn2_try, nind_current}, xnorm, ometa2);
                    if (precision_loop || *gpeucn2 <= epsgpen2) {
                        *inform = 0;
                        return;
                    }
                    if (*gpsupn <= *epsgpsn) {
                        *inform = 1;
                        return;
                    }

                    const double currprog = fprev - *f;
                    bestprog = std::max(currprog, bestprog);
                    if (currprog <= (*epsnfp) * bestprog) {
                        itnfp += 1;
                        if (itnfp >= *maxitnfp) {
                            *inform = 2;
                            return;
                        }
                    } else {
                        itnfp = 0;
                    }

                    if (*maxitngp > 0) {
                        double gpnmax = 0.0;
                        for (int i = 0; i < *maxitngp; ++i) {
                            gpnmax = std::max(gpnmax, lastgpns[i]);
                        }
                        lastgpns[*iter % (*maxitngp)] = *gpeucn2;
                        if (*gpeucn2 >= gpnmax) {
                            *inform = 3;
                            return;
                        }
                    }

                    if (*f <= *fmin) {
                        *inform = 4;
                        return;
                    }
                    if (*iter >= *maxit) {
                        *inform = 7;
                        return;
                    }
                    if (*fcnt >= *maxfc) {
                        *inform = 8;
                        return;
                    }

                    *iter += 1;
                    fprev = *f;

                    std::vector<double> x_current(n_val);
                    std::vector<double> g_current(n_val);
                    std::vector<int> ind_current(nind_current);
                    for (int i = 0; i < n_val; ++i) {
                        x_current[i] = x[i];
                        g_current[i] = g[i];
                    }
                    for (int i = 0; i < nind_current; ++i) {
                        ind_current[i] = ind[i];
                    }

                    int step_inform = -1;
                    if (gpeucn2_try > 0.0 && gieucn2_try <= ometa2 * gpeucn2_try) {
                        if (gencan_debug_enabled()) {
                            std::fprintf(stderr, "[gencan-cpp] choose=SPG iter=%d\n", *iter);
                        }
                        step_inform = run_spg_post_phase(
                            problem, options, commit_view, *iter, xnorm, sts, sty, *f, *gpeucn2,
                            x_current, g_current
                        );
                    } else {
                        if (gencan_debug_enabled()) {
                            std::fprintf(stderr, "[gencan-cpp] choose=TN iter=%d nind=%d\n", *iter, nind_current);
                        }
                        step_inform = run_tn_post_phase(
                            problem, options, commit_view, cg_schedule_seed, *iter, xnorm, sts, sty,
                            *f, *gpsupn, *gpeucn2, nind_current, x_current, g_current, ind_current
                        );
                    }

                    if (step_inform < 0) {
                        return;
                    }
                    if (step_inform != kContinueMainLoop) {
                        return;
                    }

                    xnorm = std::sqrt(norm2_kernel(n_val, x));
                    sts = 0.0;
                    sty = 0.0;
                    for (int i = 0; i < n_val; ++i) {
                        sts += s[i] * s[i];
                        sty += s[i] * y[i];
                    }

                    const ProjectedGradientSummary current_summary =
                        analyze_projected_gradient(n_val, x, g, l, u, ind);
                    gpsupn_try = current_summary.gpsupn;
                    gpeucn2_try = current_summary.gpeucn2;
                    gieucn2_try = current_summary.gieucn2;
                    nind_current = current_summary.nind;
                    *gpeucn2 = gpeucn2_try;
                    *gpsupn = gpsupn_try;
                    *inform = -1;
                }
            }

        }
    }
}
