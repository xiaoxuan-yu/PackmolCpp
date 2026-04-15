#include "gencan/types.hpp"

#include <cctype>
#include <cstdlib>
#include <string>

namespace {

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
        return GencanImplMode::kCpp;
    }
    return GencanImplMode::kCpp;
}

bool parse_env_flag_default_false(const char* name) {
    const char* env = std::getenv(name);
    if (env == nullptr) {
        return false;
    }

    std::string value(env);
    for (char& ch : value) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value == "1" || value == "true" || value == "on" || value == "yes";
}

} // namespace

GencanImplMode active_impl_mode() {
    static const GencanImplMode mode = parse_impl_mode();
    return mode;
}

bool gencan_debug_enabled() {
    return parse_env_flag_default_false("PACKMOL_GENCAN_DEBUG");
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

bool use_cpp_eval_kernel() {
    return parse_env_flag_default_false("PACKMOL_GENCAN_EVAL_CPP");
}
