#include "gencan/api.hpp"

#include <cstdint>

extern "C" void packmol_eval_meta_state_c(
    int* maxrest,
    int* mrperatom,
    int* ntype,
    int* ntotmol,
    int* ntotat,
    int* nfixedat,
    int* ntype_with_fixed,
    int* lcellfirst,
    int* ncells,
    double* scale,
    double* scale2,
    double* pbc_min,
    double* pbc_length,
    double* cell_length,
    bool* init1,
    bool* move,
    bool* fix,
    bool* using_pbc
);

extern "C" void packmol_eval_numeric_pointers_c(
    void** nmols_ptr,
    void** natoms_ptr,
    void** idfirst_ptr,
    void** nratom_ptr,
    void** iratom_ptr,
    void** ityperest_ptr,
    void** ibmol_ptr,
    void** ibtype_ptr,
    void** latomnext_ptr,
    void** latomfirst_ptr,
    void** lcellnext_ptr,
    void** xcart_ptr,
    void** coor_ptr,
    void** restpars_ptr,
    void** radius_ptr,
    void** radius_ini_ptr,
    void** fscale_ptr,
    void** short_radius_ptr,
    void** short_radius_scale_ptr,
    void** gxcar_ptr,
    void** fdist_atom_ptr,
    void** frest_atom_ptr
);

extern "C" void packmol_eval_copy_type_flags_c(int* comptype_out);
extern "C" void packmol_eval_copy_atom_flags_c(int* fixedatom_out, int* use_short_radius_out);

bool load_packmol_eval_state_cpp(PackmolEvalStateSnapshot* state) {
    if (state == nullptr) {
        return false;
    }

    PackmolEvalMetaState meta{};
    packmol_eval_meta_state_c(
        &meta.maxrest,
        &meta.mrperatom,
        &meta.ntype,
        &meta.ntotmol,
        &meta.ntotat,
        &meta.nfixedat,
        &meta.ntype_with_fixed,
        &meta.lcellfirst,
        meta.ncells,
        &meta.scale,
        &meta.scale2,
        meta.pbc_min,
        meta.pbc_length,
        meta.cell_length,
        &meta.init1,
        &meta.move,
        &meta.fix,
        &meta.using_pbc
    );

    PackmolEvalRawPointers raw{};
    void* nmols_ptr = nullptr;
    void* natoms_ptr = nullptr;
    void* idfirst_ptr = nullptr;
    void* nratom_ptr = nullptr;
    void* iratom_ptr = nullptr;
    void* ityperest_ptr = nullptr;
    void* ibmol_ptr = nullptr;
    void* ibtype_ptr = nullptr;
    void* latomnext_ptr = nullptr;
    void* latomfirst_ptr = nullptr;
    void* lcellnext_ptr = nullptr;
    void* xcart_ptr = nullptr;
    void* coor_ptr = nullptr;
    void* restpars_ptr = nullptr;
    void* radius_ptr = nullptr;
    void* radius_ini_ptr = nullptr;
    void* fscale_ptr = nullptr;
    void* short_radius_ptr = nullptr;
    void* short_radius_scale_ptr = nullptr;
    void* gxcar_ptr = nullptr;
    void* fdist_atom_ptr = nullptr;
    void* frest_atom_ptr = nullptr;
    packmol_eval_numeric_pointers_c(
        &nmols_ptr,
        &natoms_ptr,
        &idfirst_ptr,
        &nratom_ptr,
        &iratom_ptr,
        &ityperest_ptr,
        &ibmol_ptr,
        &ibtype_ptr,
        &latomnext_ptr,
        &latomfirst_ptr,
        &lcellnext_ptr,
        &xcart_ptr,
        &coor_ptr,
        &restpars_ptr,
        &radius_ptr,
        &radius_ini_ptr,
        &fscale_ptr,
        &short_radius_ptr,
        &short_radius_scale_ptr,
        &gxcar_ptr,
        &fdist_atom_ptr,
        &frest_atom_ptr
    );
    raw.nmols = static_cast<const int*>(nmols_ptr);
    raw.natoms = static_cast<const int*>(natoms_ptr);
    raw.idfirst = static_cast<const int*>(idfirst_ptr);
    raw.nratom = static_cast<const int*>(nratom_ptr);
    raw.iratom = static_cast<const int*>(iratom_ptr);
    raw.ityperest = static_cast<const int*>(ityperest_ptr);
    raw.ibmol = static_cast<int*>(ibmol_ptr);
    raw.ibtype = static_cast<int*>(ibtype_ptr);
    raw.latomnext = static_cast<int*>(latomnext_ptr);
    raw.latomfirst = static_cast<int*>(latomfirst_ptr);
    raw.lcellnext = static_cast<int*>(lcellnext_ptr);
    raw.xcart = static_cast<double*>(xcart_ptr);
    raw.coor = static_cast<const double*>(coor_ptr);
    raw.restpars = static_cast<const double*>(restpars_ptr);
    raw.radius = static_cast<const double*>(radius_ptr);
    raw.radius_ini = static_cast<const double*>(radius_ini_ptr);
    raw.fscale = static_cast<const double*>(fscale_ptr);
    raw.short_radius = static_cast<const double*>(short_radius_ptr);
    raw.short_radius_scale = static_cast<const double*>(short_radius_scale_ptr);
    raw.gxcar = static_cast<double*>(gxcar_ptr);
    raw.fdist_atom = static_cast<double*>(fdist_atom_ptr);
    raw.frest_atom = static_cast<double*>(frest_atom_ptr);

    state->meta = meta;
    state->raw = raw;
    state->comptype_flags.assign(meta.ntype, 0);
    state->fixedatom_flags.assign(meta.ntotat, 0);
    state->use_short_radius_flags.assign(meta.ntotat, 0);

    if (meta.ntype > 0) {
        packmol_eval_copy_type_flags_c(state->comptype_flags.data());
    }
    if (meta.ntotat > 0) {
        packmol_eval_copy_atom_flags_c(
            state->fixedatom_flags.data(), state->use_short_radius_flags.data()
        );
    }

    return raw.nmols != nullptr && raw.natoms != nullptr && raw.idfirst != nullptr &&
           raw.nratom != nullptr && raw.iratom != nullptr && raw.ityperest != nullptr &&
           raw.xcart != nullptr && raw.coor != nullptr && raw.restpars != nullptr &&
           raw.radius != nullptr && raw.radius_ini != nullptr && raw.fscale != nullptr &&
           raw.short_radius != nullptr && raw.short_radius_scale != nullptr &&
           raw.gxcar != nullptr;
}
