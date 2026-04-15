#include "gencan/api.hpp"

#include <algorithm>
#include <cmath>

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

namespace {

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

} // namespace

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

void run_tnls_context_cpp(const GencanTnlsContext& context) {
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_tnls_fortran_c(
                context.nind, context.ind, context.n, context.x, context.m, context.lambda,
                context.rho, context.l, context.u, context.f, context.g, context.d, context.amax,
                context.rbdtype, context.rbdind, context.nint, context.next, context.mininterp,
                context.maxextrap, context.fmin, context.maxfc, context.gtype, context.iprint,
                context.fcnt, context.gcnt, context.intcnt, context.exgcnt, context.exbcnt,
                context.inform, context.xplus, context.xtmp, context.xbext, context.gamma,
                context.beta, context.sigma1, context.sigma2, context.sterel, context.steabs,
                context.epsrel, context.epsabs, context.infrel, context.infabs
            );
            return;
        case GencanImplMode::kCpp:
            if (tnls_cpp_subset(
                    context.nind, context.ind, context.n, context.x, context.m, context.lambda,
                    context.rho, context.l, context.u, context.f, context.g, context.d, context.amax,
                    context.rbdtype, context.rbdind, context.nint, context.next, context.mininterp,
                    context.maxextrap, context.fmin, context.maxfc, context.gtype, context.fcnt,
                    context.gcnt, context.intcnt, context.exgcnt, context.exbcnt, context.inform,
                    context.xplus, context.xtmp, context.xbext, context.gamma, context.beta,
                    context.sigma1, context.sigma2, context.sterel, context.steabs, context.epsrel,
                    context.epsabs)) {
                return;
            }
            return;
    }
}

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
) {
    run_tnls_context_cpp(GencanTnlsContext{
        nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, rbdind, nint, next,
        mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, gcnt, intcnt, exgcnt, exbcnt,
        inform, xplus, xtmp, xbext, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel,
        epsabs, infrel, infabs
    });
}
