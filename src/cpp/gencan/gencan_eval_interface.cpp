#include "gencan/api.hpp"

#include <algorithm>
#include <cmath>

namespace {

thread_local int g_init1_phase_active = 0;

} // namespace

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
    (void)g_init1_phase_active;
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
    *flag = 0;
    for (int j = 0; j < *n; ++j) {
        const double xj = x[j];
        const double step = std::max(*steabs, (*sterel) * std::abs(xj));

        x[j] = xj + step;
        double fplus = 0.0;
        packmol_evalal_fortran_c(n, x, m, lambda, rho, &fplus, flag);
        if (*flag < 0) {
            x[j] = xj;
            return;
        }

        x[j] = xj - step;
        double fminus = 0.0;
        packmol_evalal_fortran_c(n, x, m, lambda, rho, &fminus, flag);
        if (*flag < 0) {
            x[j] = xj;
            return;
        }

        g[j] = (fplus - fminus) / (2.0 * step);
        x[j] = xj;
    }
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

extern "C" void packmol_evalal_cpp_direct_c(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* flag,
    int* used_cpp
) {
    const bool ok = eval_objective_full_cpp_direct(n, x, m, lambda, rho, f, flag);
    if (used_cpp != nullptr) {
        *used_cpp = ok ? 1 : 0;
    }
    if (!ok && flag != nullptr) {
        *flag = -1;
    }
}
