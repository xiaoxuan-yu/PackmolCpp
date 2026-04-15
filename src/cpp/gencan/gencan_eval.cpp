#include "gencan/api.hpp"

#include <algorithm>
#include <cmath>

extern "C" void packmol_evalal_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
);

extern "C" void packmol_evalnal_fortran_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    int* inform
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
    int* inform
);

extern "C" void packmol_evalhd_fortran_c(const int* n);

extern "C" void packmol_precision_state_fortran_c(
    double* precision_value,
    double* fdist_value,
    double* frest_value
);

void eval_objective_full_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
) {
    if (try_eval_objective_full_cpp(n, x, m, lambda, rho, f, inform)) {
        return;
    }
    packmol_evalal_fortran_c(n, x, m, lambda, rho, f, inform);
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
    eval_objective_full_cpp(n, x, m, lambda, rho, f, inform);
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
    eval_gradient_full_cpp(n, x, m, lambda, rho, gtype, g, sterel, steabs, inform);
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
        return;
    }

    *inform = 0;
    for (int j = 0; j < *n; ++j) {
        const double xj = x[j];
        const double step = std::max(*steabs, (*sterel) * std::abs(xj));

        x[j] = xj + step;
        double fplus = 0.0;
        eval_objective_full_cpp(n, x, m, lambda, rho, &fplus, inform);
        if (*inform < 0) {
            x[j] = xj;
            return;
        }

        x[j] = xj - step;
        double fminus = 0.0;
        eval_objective_full_cpp(n, x, m, lambda, rho, &fminus, inform);
        if (*inform < 0) {
            x[j] = xj;
            return;
        }

        g[j] = (fplus - fminus) / (2.0 * step);
        x[j] = xj;
    }
}

PackmolPrecisionState read_packmol_precision_state_cpp() {
    PackmolPrecisionState state{0.0, 0.0, 0.0};
    packmol_precision_state_fortran_c(
        &state.precision_value, &state.fdist_value, &state.frest_value
    );
    return state;
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
    eval_objective_full_cpp(n, x, m, lambda, rho, &f_unused, eval_flag);
    if (*eval_flag < 0) {
        return false;
    }
    const PackmolPrecisionState state = read_packmol_precision_state_cpp();
    return state.fdist_value < state.precision_value && state.frest_value < state.precision_value;
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

    g_evalhd_fortran_call_count += 1;
    packmol_evalhd_fortran_c(n);

    shrink_inplace(nind_val, ind, x);
    shrink_inplace(nind_val, ind, d);
    shrink_inplace(nind_val, ind, g);
    shrink_inplace(nind_val, ind, hd);
}
