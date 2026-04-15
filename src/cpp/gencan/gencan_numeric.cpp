#include "gencan/api.hpp"

#include <cmath>

namespace {

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

} // namespace

double norm2_kernel(const int n, const double* x) {
    if (use_cpp_numeric_kernel()) {
        return norm2sq_cpp_stable(n, x);
    }
    const double nrm = hsldnrm2_cpp_legacy(n, x);
    return nrm * nrm;
}
