#include "gencan/api.hpp"

#include <algorithm>
#include <cmath>

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

int finalize_post_step_inform_cpp(
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
