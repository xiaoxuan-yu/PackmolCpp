#include "gencan/api.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

thread_local int g_cg_active_call_id = 0;

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

double dot_cpp_legacy(const int n, const double* a, const double* b) {
    double val = 0.0;
    for (int i = 0; i < n; ++i) {
        val += a[i] * b[i];
    }
    return val;
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

double norm2_kernel_local(const int n, const double* x) {
    if (use_cpp_numeric_kernel()) {
        return norm2sq_cpp_stable(n, x);
    }
    const double nrm = hsldnrm2_cpp_legacy(n, x);
    return nrm * nrm;
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
    (void)ind;
    (void)n;
    (void)x;
    (void)m;
    (void)lambda;
    (void)rho;
    (void)iprint;
    const bool debug = gencan_debug_enabled();
    *rbdtype = 0;
    *rbdind = 0;
    auto norm2 = [&](const double* v) { return norm2_kernel_local(nind_val, v); };
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
                    return;
                }
            }
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
            return;
        }

        if (alpha == amax1 || alpha == amax1n) {
            *rbdtype = 0;
            *rbdind = 0;
            *inform = 1;
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
                return;
            }
        } else {
            itnqmp = 0;
        }

        (void)debug;
    }
}

} // namespace

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

void run_cg_context_cpp(const GencanCgContext& context) {
    static int cg_call_counter = 0;
    cg_call_counter += 1;
    switch (active_impl_mode()) {
        case GencanImplMode::kFortran:
            packmol_cg_fortran_c(
                context.nind, context.ind, context.n, context.x, context.m, context.lambda,
                context.rho, context.g, context.delta, context.l, context.u, context.eps,
                context.epsnqmp, context.maxitnqmp, context.maxit, context.nearlyq, context.gtype,
                context.htvtype, context.trtype, context.iprint, context.ncomp, context.s,
                context.iter, context.rbdtype, context.rbdind, context.inform, context.w,
                context.y, context.r, context.d, context.sprev, context.theta, context.sterel,
                context.steabs, context.epsrel, context.epsabs, context.infrel, context.infabs
            );
            return;
        case GencanImplMode::kCpp:
            g_cg_active_call_id = cg_call_counter;
            cg_cpp_full(
                context.nind, context.ind, context.n, context.x, context.m, context.lambda,
                context.rho, context.g, context.delta, context.l, context.u, context.eps,
                context.epsnqmp, context.maxitnqmp, context.maxit, context.nearlyq, context.gtype,
                context.htvtype, context.trtype, context.iprint, context.theta, context.sterel,
                context.steabs, context.epsrel, context.epsabs, context.infabs, context.s,
                context.w, context.y, context.r, context.d, context.sprev, context.iter,
                context.rbdtype, context.rbdind, context.inform
            );
            g_cg_active_call_id = 0;
            return;
    }
}

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
) {
    run_cg_context_cpp(GencanCgContext{
        nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp, maxitnqmp, maxit,
        nearlyq, gtype, htvtype, trtype, iprint, ncomp, s, iter, rbdtype, rbdind, inform, w,
        y, r, d, sprev, theta, sterel, steabs, epsrel, epsabs, infrel, infabs
    });
}
