#include "gencan/api.hpp"

#include <algorithm>
#include <cmath>

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

} // namespace

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
    eval_objective_full_cpp(n, xtrial, m, lambda, rho, &ftrial, inform);
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

        eval_objective_full_cpp(&n_val, xtrial, &m_val, lambda, rho, &ftrial, inform);
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

void run_spgls_context_cpp(const GencanSpglsContext& context) {
    const bool use_fortran = active_impl_mode() == GencanImplMode::kFortran || !use_cpp_spgls_kernel();
    if (use_fortran) {
        packmol_spgls_fortran_c(
            context.n, context.x, context.m, context.lambda, context.rho, context.f, context.g,
            context.l, context.u, context.lamspg, context.nint, context.mininterp, context.fmin,
            context.maxfc, context.iprint, context.fcnt, context.inform, context.xtrial, context.d,
            context.gamma, context.sigma1, context.sigma2, context.sterel, context.steabs,
            context.epsrel, context.epsabs, context.infrel, context.infabs
        );
        return;
    }
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_spgls_fortran_c(
                context.n, context.x, context.m, context.lambda, context.rho, context.f, context.g,
                context.l, context.u, context.lamspg, context.nint, context.mininterp, context.fmin,
                context.maxfc, context.iprint, context.fcnt, context.inform, context.xtrial, context.d,
                context.gamma, context.sigma1, context.sigma2, context.sterel, context.steabs,
                context.epsrel, context.epsabs, context.infrel, context.infabs
            );
            return;
        case GencanImplMode::kCpp:
            spgls_cpp(
                context.n, context.x, context.m, context.lambda, context.rho, context.f, context.g,
                context.l, context.u, context.lamspg, context.nint, context.mininterp, context.fmin,
                context.maxfc, context.fcnt, context.inform, context.xtrial, context.d, context.gamma,
                context.sigma1, context.sigma2, context.epsrel, context.epsabs
            );
            return;
    }
}

void packmol_gencan_spgls_bridge_impl(
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
    run_spgls_context_cpp(GencanSpglsContext{
        n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, fmin, maxfc, iprint,
        fcnt, inform, xtrial, d, gamma, sigma1, sigma2, sterel, steabs, epsrel, epsabs,
        infrel, infabs
    });
}
