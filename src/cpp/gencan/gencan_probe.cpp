#include "gencan/api.hpp"

int packmol_gencan_cpp_probe_impl(const int n, const double* x, double* fx_out) {
    if (n < 0 || x == nullptr || fx_out == nullptr) {
        return -1;
    }

    *fx_out = 0.0;
    return 0;
}

int packmol_gencan_impl_mode_cpp() {
    return static_cast<int>(active_impl_mode());
}

void packmol_gencan_reset_evalhd_call_count_cpp() {
    g_evalhd_fortran_call_count = 0;
}

int packmol_gencan_evalhd_call_count_cpp() {
    return g_evalhd_fortran_call_count;
}
