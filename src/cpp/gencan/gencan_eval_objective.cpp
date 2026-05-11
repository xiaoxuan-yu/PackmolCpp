#include "gencan/api.hpp"

#include <algorithm>
#include <array>
#include <cmath>

extern "C" double packmol_fdist;
extern "C" double packmol_frest;
extern "C" double packmol_pair_penalty_sum;
extern "C" double packmol_constraint_penalty_sum;
extern "C" int packmol_pair_penalty_count;
extern "C" int packmol_constraint_penalty_count;
extern "C" void packmol_eval_set_lcellfirst_c(const int* lcellfirst);

namespace {

int latomfirst_index(const PackmolEvalStateSnapshot& state, const int i, const int j, const int k) {
    return (i - 1) + (j - 1) * state.meta.ncells[0] +
           (k - 1) * state.meta.ncells[0] * state.meta.ncells[1];
}

double& xcart_ref(const PackmolEvalStateSnapshot& state, const int icart, const int axis) {
    return state.raw.xcart[(icart - 1) + (axis - 1) * state.meta.ntotat];
}

double xcart_value(const PackmolEvalStateSnapshot& state, const int icart, const int axis) {
    return state.raw.xcart[(icart - 1) + (axis - 1) * state.meta.ntotat];
}

double coor_value(const PackmolEvalStateSnapshot& state, const int iatom, const int axis) {
    return state.raw.coor[(iatom - 1) + (axis - 1) * state.meta.ntotat];
}

double restpars_value(const PackmolEvalStateSnapshot& state, const int irest, const int axis) {
    return state.raw.restpars[(irest - 1) + (axis - 1) * state.meta.maxrest];
}

int iratom_value(const PackmolEvalStateSnapshot& state, const int icart, const int iratcount) {
    return state.raw.iratom[(icart - 1) + (iratcount - 1) * state.meta.ntotat];
}

int cell_ind_cpp(const int icell, const int ncells) {
    return ((icell - 1 + ncells) % ncells) + 1;
}

double v_in_box_cpp(const double v, const double pbc_min, const double pbc_length) {
    return pbc_min + std::fmod(std::fmod(v - pbc_min, pbc_length) + pbc_length, pbc_length);
}

double delta_vector_cpp(const double v1, const double v2, const double pbc_length) {
    double delta = v1 - v2;
    delta -= pbc_length * std::nearbyint(delta / pbc_length);
    return delta;
}

void setcell_cpp(const PackmolEvalStateSnapshot& state, const double* pos, int* cell) {
    for (int axis = 0; axis < 3; ++axis) {
        const double xt = v_in_box_cpp(
            pos[axis], state.meta.pbc_min[axis], state.meta.pbc_length[axis]
        );
        int c = static_cast<int>(
                    std::floor((xt - state.meta.pbc_min[axis]) / state.meta.cell_length[axis])) +
                1;
        cell[axis] = std::min(c, state.meta.ncells[axis]);
    }
}

int index_cell_cpp(const PackmolEvalStateSnapshot& state, const int* cell) {
    return (cell[0] - 1) * state.meta.ncells[1] * state.meta.ncells[2] +
           (cell[1] - 1) * state.meta.ncells[2] + (cell[2] - 1) + 1;
}

void icell_to_cell_cpp(const PackmolEvalStateSnapshot& state, const int icell, int* cell) {
    int iicell = icell;
    cell[2] = iicell % state.meta.ncells[2];
    if (cell[2] == 0) {
        cell[2] = state.meta.ncells[2];
    }
    iicell = (iicell - cell[2]) / state.meta.ncells[2] + 1;
    cell[1] = iicell % state.meta.ncells[1];
    if (cell[1] == 0) {
        cell[1] = state.meta.ncells[1];
    }
    iicell = (iicell - cell[1]) / state.meta.ncells[1] + 1;
    cell[0] = iicell % state.meta.ncells[0];
    if (cell[0] == 0) {
        cell[0] = state.meta.ncells[0];
    }
}

void resetcells_cpp(const PackmolEvalStateSnapshot& state, int* lcellfirst) {
    const int ncells_total =
        state.meta.ncells[0] * state.meta.ncells[1] * state.meta.ncells[2];
    for (int i = 0; i < ncells_total; ++i) {
        state.raw.lcellnext[i] = 0;
        state.raw.latomfirst[i] = 0;
    }
    for (int i = 0; i < state.meta.ntotat; ++i) {
        state.raw.latomnext[i] = 0;
    }
    *lcellfirst = 0;

    if (!state.meta.fix) {
        return;
    }

    int icart = state.meta.ntotat - state.meta.nfixedat;
    for (int iftype = state.meta.ntype + 1; iftype <= state.meta.ntype_with_fixed; ++iftype) {
        for (int ifatom = 1; ifatom <= state.raw.natoms[iftype - 1]; ++ifatom) {
            icart += 1;
            const double pos[3] = {
                xcart_value(state, icart, 1),
                xcart_value(state, icart, 2),
                xcart_value(state, icart, 3),
            };
            int cell[3] = {0, 0, 0};
            setcell_cpp(state, pos, cell);
            const int latomfirst_idx = latomfirst_index(state, cell[0], cell[1], cell[2]);
            const bool empty_cell = state.raw.latomfirst[latomfirst_idx] == 0;
            state.raw.latomnext[icart - 1] = state.raw.latomfirst[latomfirst_idx];
            state.raw.latomfirst[latomfirst_idx] = icart;
            if (empty_cell) {
                const int icell = index_cell_cpp(state, cell);
                state.raw.lcellnext[icell - 1] = *lcellfirst;
                *lcellfirst = icell;
            }
        }
    }
}

void eulerrmat_cpp(
    const double beta,
    const double gama,
    const double teta,
    std::array<double, 3>& v1,
    std::array<double, 3>& v2,
    std::array<double, 3>& v3
) {
    const double cb = std::cos(beta);
    const double sb = std::sin(beta);
    const double cg = std::cos(gama);
    const double sg = std::sin(gama);
    const double ct = std::cos(teta);
    const double st = std::sin(teta);

    v1[0] = -sb * sg * ct + cb * cg;
    v1[1] = -sb * cg * ct - cb * sg;
    v1[2] = sb * st;

    v2[0] = cb * sg * ct + sb * cg;
    v2[1] = cb * cg * ct - sb * sg;
    v2[2] = -cb * st;

    v3[0] = sg * st;
    v3[1] = cg * st;
    v3[2] = ct;
}

void compcart_cpp(
    double* out,
    const std::array<double, 3>& xcm,
    const std::array<double, 3>& xref,
    const std::array<double, 3>& v1,
    const std::array<double, 3>& v2,
    const std::array<double, 3>& v3
) {
    out[0] = xcm[0] + xref[0] * v1[0] + xref[1] * v2[0] + xref[2] * v3[0];
    out[1] = xcm[1] + xref[0] * v1[1] + xref[1] * v2[1] + xref[2] * v3[1];
    out[2] = xcm[2] + xref[0] * v1[2] + xref[1] * v2[2] + xref[2] * v3[2];
}

double comprest_atom_cpp(const PackmolEvalStateSnapshot& state, const int icart) {
    const double scale = state.meta.scale;
    const double scale2 = state.meta.scale2;
    const double x = xcart_value(state, icart, 1);
    const double y = xcart_value(state, icart, 2);
    const double z = xcart_value(state, icart, 3);

    double f = 0.0;
    for (int iratcount = 1; iratcount <= state.raw.nratom[icart - 1]; ++iratcount) {
        const int irest = iratom_value(state, icart, iratcount);
        const int type = state.raw.ityperest[irest - 1];
        if (type == 2) {
            const double clength = restpars_value(state, irest, 4);
            const double xmin = restpars_value(state, irest, 1);
            const double ymin = restpars_value(state, irest, 2);
            const double zmin = restpars_value(state, irest, 3);
            const double xmax = xmin + clength;
            const double ymax = ymin + clength;
            const double zmax = zmin + clength;
            double a1 = std::min(x - xmin, 0.0);
            double a2 = std::min(y - ymin, 0.0);
            double a3 = std::min(z - zmin, 0.0);
            f += scale * (a1 * a1 + a2 * a2 + a3 * a3);
            a1 = std::max(x - xmax, 0.0);
            a2 = std::max(y - ymax, 0.0);
            a3 = std::max(z - zmax, 0.0);
            f += scale * (a1 * a1 + a2 * a2 + a3 * a3);
        } else if (type == 3) {
            const double xmin = restpars_value(state, irest, 1);
            const double ymin = restpars_value(state, irest, 2);
            const double zmin = restpars_value(state, irest, 3);
            const double xmax = restpars_value(state, irest, 4);
            const double ymax = restpars_value(state, irest, 5);
            const double zmax = restpars_value(state, irest, 6);
            double a1 = std::min(x - xmin, 0.0);
            double a2 = std::min(y - ymin, 0.0);
            double a3 = std::min(z - zmin, 0.0);
            f += scale * (a1 * a1 + a2 * a2 + a3 * a3);
            a1 = std::max(x - xmax, 0.0);
            a2 = std::max(y - ymax, 0.0);
            a3 = std::max(z - zmax, 0.0);
            f += scale * (a1 * a1 + a2 * a2 + a3 * a3);
        } else if (type == 4) {
            const double w =
                (x - restpars_value(state, irest, 1)) * (x - restpars_value(state, irest, 1)) +
                (y - restpars_value(state, irest, 2)) * (y - restpars_value(state, irest, 2)) +
                (z - restpars_value(state, irest, 3)) * (z - restpars_value(state, irest, 3)) -
                restpars_value(state, irest, 4) * restpars_value(state, irest, 4);
            const double a1 = std::max(w, 0.0);
            f += scale2 * a1 * a1;
        } else if (type == 5) {
            const double a1 =
                (x - restpars_value(state, irest, 1)) * (x - restpars_value(state, irest, 1)) /
                std::pow(restpars_value(state, irest, 4), 2);
            const double a2 =
                (y - restpars_value(state, irest, 2)) * (y - restpars_value(state, irest, 2)) /
                std::pow(restpars_value(state, irest, 5), 2);
            const double a3 =
                (z - restpars_value(state, irest, 3)) * (z - restpars_value(state, irest, 3)) /
                std::pow(restpars_value(state, irest, 6), 2);
            const double a4 = std::pow(restpars_value(state, irest, 7), 2);
            const double w = a1 + a2 + a3 - a4;
            const double p = std::max(w, 0.0);
            f += scale2 * p * p;
        } else if (type == 6 || type == 7) {
            const double xmin = restpars_value(state, irest, 1);
            const double ymin = restpars_value(state, irest, 2);
            const double zmin = restpars_value(state, irest, 3);
            const double xmax = (type == 6) ? xmin + restpars_value(state, irest, 4)
                                            : restpars_value(state, irest, 4);
            const double ymax = (type == 6) ? ymin + restpars_value(state, irest, 4)
                                            : restpars_value(state, irest, 5);
            const double zmax = (type == 6) ? zmin + restpars_value(state, irest, 4)
                                            : restpars_value(state, irest, 6);
            double a1 = 0.0;
            double a2 = 0.0;
            double a3 = 0.0;
            if ((x > xmin && x < xmax) && (y > ymin && y < ymax) && (z > zmin && z < zmax)) {
                const double xmed = (xmax - xmin) / 2.0;
                const double ymed = (ymax - ymin) / 2.0;
                const double zmed = (zmax - zmin) / 2.0;
                a1 = (x <= xmed) ? (x - xmin) : (xmax - x);
                a2 = (y <= ymed) ? (y - ymin) : (ymax - y);
                a3 = (z <= zmed) ? (z - zmin) : (zmax - z);
            }
            f += scale * (a1 + a2 + a3);
        } else if (type == 8) {
            const double w =
                (x - restpars_value(state, irest, 1)) * (x - restpars_value(state, irest, 1)) +
                (y - restpars_value(state, irest, 2)) * (y - restpars_value(state, irest, 2)) +
                (z - restpars_value(state, irest, 3)) * (z - restpars_value(state, irest, 3)) -
                restpars_value(state, irest, 4) * restpars_value(state, irest, 4);
            const double a1 = std::min(w, 0.0);
            f += scale2 * a1 * a1;
        } else if (type == 9) {
            const double a1 =
                (x - restpars_value(state, irest, 1)) * (x - restpars_value(state, irest, 1)) /
                std::pow(restpars_value(state, irest, 4), 2);
            const double a2 =
                (y - restpars_value(state, irest, 2)) * (y - restpars_value(state, irest, 2)) /
                std::pow(restpars_value(state, irest, 5), 2);
            const double a3 =
                (z - restpars_value(state, irest, 3)) * (z - restpars_value(state, irest, 3)) /
                std::pow(restpars_value(state, irest, 6), 2);
            const double a4 = std::pow(restpars_value(state, irest, 7), 2);
            const double w = a1 + a2 + a3 - a4;
            const double p = std::min(w, 0.0);
            f += p * p;
        } else if (type == 10 || type == 11) {
            const double w = restpars_value(state, irest, 1) * x + restpars_value(state, irest, 2) * y +
                             restpars_value(state, irest, 3) * z - restpars_value(state, irest, 4);
            const double p = (type == 10) ? std::min(w, 0.0) : std::max(w, 0.0);
            f += scale * p * p;
        } else if (type == 12 || type == 13) {
            const double a1 = x - restpars_value(state, irest, 1);
            const double a2 = y - restpars_value(state, irest, 2);
            const double a3 = z - restpars_value(state, irest, 3);
            const double vnorm = std::sqrt(
                std::pow(restpars_value(state, irest, 4), 2) +
                std::pow(restpars_value(state, irest, 5), 2) +
                std::pow(restpars_value(state, irest, 6), 2)
            );
            const double v1 = restpars_value(state, irest, 4) / vnorm;
            const double v2 = restpars_value(state, irest, 5) / vnorm;
            const double v3 = restpars_value(state, irest, 6) / vnorm;
            const double w = v1 * a1 + v2 * a2 + v3 * a3;
            const double d = std::pow(a1 - v1 * w, 2) + std::pow(a2 - v2 * w, 2) +
                             std::pow(a3 - v3 * w, 2);
            if (type == 12) {
                f += scale2 * (std::pow(std::max(-w, 0.0), 2) +
                               std::pow(std::max(w - restpars_value(state, irest, 9), 0.0), 2) +
                               std::pow(std::max(d - std::pow(restpars_value(state, irest, 7), 2), 0.0), 2));
            } else {
                f += scale2 * (std::pow(std::min(-w, 0.0), 2) *
                               std::pow(std::min(w - restpars_value(state, irest, 9), 0.0), 2) *
                               std::pow(std::min(d - std::pow(restpars_value(state, irest, 7), 2), 0.0), 2));
            }
        } else if (type == 14 || type == 15) {
            const double a1 = -std::pow(x - restpars_value(state, irest, 1), 2) /
                              (2.0 * std::pow(restpars_value(state, irest, 3), 2));
            const double a2 = -std::pow(y - restpars_value(state, irest, 2), 2) /
                              (2.0 * std::pow(restpars_value(state, irest, 4), 2));
            double w = -(z - restpars_value(state, irest, 5));
            if (a1 + a2 > -50.0) {
                w = restpars_value(state, irest, 6) * std::exp(a1 + a2) -
                    (z - restpars_value(state, irest, 5));
            }
            const double p = (type == 14) ? std::max(w, 0.0) : std::min(w, 0.0);
            f += scale * p * p;
        }
    }
    return f;
}

double fparc_cpp(const PackmolEvalStateSnapshot& state, const int icart, int firstjcart) {
    double fparc = 0.0;
    int jcart = firstjcart;
    while (jcart > 0) {
        if (state.comptype_flags[state.raw.ibtype[jcart - 1] - 1] == 0) {
            jcart = state.raw.latomnext[jcart - 1];
            continue;
        }
        if (state.raw.ibmol[icart - 1] == state.raw.ibmol[jcart - 1] &&
            state.raw.ibtype[icart - 1] == state.raw.ibtype[jcart - 1]) {
            jcart = state.raw.latomnext[jcart - 1];
            continue;
        }
        if (state.fixedatom_flags[icart - 1] != 0 && state.fixedatom_flags[jcart - 1] != 0) {
            jcart = state.raw.latomnext[jcart - 1];
            continue;
        }

        const double dx = delta_vector_cpp(
            xcart_value(state, icart, 1), xcart_value(state, jcart, 1), state.meta.pbc_length[0]
        );
        const double dy = delta_vector_cpp(
            xcart_value(state, icart, 2), xcart_value(state, jcart, 2), state.meta.pbc_length[1]
        );
        const double dz = delta_vector_cpp(
            xcart_value(state, icart, 3), xcart_value(state, jcart, 3), state.meta.pbc_length[2]
        );
        const double datom = dx * dx + dy * dy + dz * dz;

        double pair_penalty = 0.0;
        double tol = std::pow(state.raw.radius[icart - 1] + state.raw.radius[jcart - 1], 2);
        if (datom < tol) {
            pair_penalty += state.raw.fscale[icart - 1] * state.raw.fscale[jcart - 1] *
                            std::pow(datom - tol, 2);
            if (state.use_short_radius_flags[icart - 1] != 0 ||
                state.use_short_radius_flags[jcart - 1] != 0) {
                const double short_tol = std::pow(
                    state.raw.short_radius[icart - 1] + state.raw.short_radius[jcart - 1], 2
                );
                if (datom < short_tol) {
                    const double short_tol_penalty = datom - short_tol;
                    double short_tol_scale = std::sqrt(
                        state.raw.short_radius_scale[icart - 1] *
                        state.raw.short_radius_scale[jcart - 1]
                    );
                    short_tol_scale *= (tol * tol) / (short_tol * short_tol);
                    pair_penalty += state.raw.fscale[icart - 1] * state.raw.fscale[jcart - 1] *
                                    short_tol_scale * short_tol_penalty * short_tol_penalty;
                }
            }
            fparc += pair_penalty;
            packmol_pair_penalty_sum += pair_penalty;
            packmol_pair_penalty_count += 1;
        }

        tol = std::pow(state.raw.radius_ini[icart - 1] + state.raw.radius_ini[jcart - 1], 2);
        packmol_fdist = std::max(packmol_fdist, tol - datom);
        if (state.meta.move && state.raw.fdist_atom != nullptr) {
            state.raw.fdist_atom[icart - 1] =
                std::max(state.raw.fdist_atom[icart - 1], tol - datom);
            state.raw.fdist_atom[jcart - 1] =
                std::max(state.raw.fdist_atom[jcart - 1], tol - datom);
        }

        jcart = state.raw.latomnext[jcart - 1];
    }
    return fparc;
}

bool eval_objective_full_cpp_kernel(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
) {
    (void)m;
    (void)lambda;
    (void)rho;
    if (n == nullptr || x == nullptr || f == nullptr || inform == nullptr) {
        return false;
    }

    PackmolEvalStateSnapshot state;
    if (!load_packmol_eval_state_cpp(&state)) {
        return false;
    }
    if (*n != state.meta.ntotmol * 6) {
        return false;
    }

    *inform = 0;
    *f = 0.0;
    packmol_fdist = 0.0;
    packmol_frest = 0.0;
    packmol_pair_penalty_sum = 0.0;
    packmol_constraint_penalty_sum = 0.0;
    packmol_pair_penalty_count = 0;
    packmol_constraint_penalty_count = 0;

    int lcellfirst = state.meta.lcellfirst;
    if (!state.meta.init1) {
        resetcells_cpp(state, &lcellfirst);
    }

    int ilubar = 0;
    int ilugan = state.meta.ntotmol * 3;
    int icart = 0;

    for (int itype = 1; itype <= state.meta.ntype; ++itype) {
        if (state.comptype_flags[itype - 1] == 0) {
            icart += state.raw.nmols[itype - 1] * state.raw.natoms[itype - 1];
            continue;
        }

        for (int imol = 1; imol <= state.raw.nmols[itype - 1]; ++imol) {
            const std::array<double, 3> xcm = {x[ilubar], x[ilubar + 1], x[ilubar + 2]};
            std::array<double, 3> v1{};
            std::array<double, 3> v2{};
            std::array<double, 3> v3{};
            eulerrmat_cpp(x[ilugan], x[ilugan + 1], x[ilugan + 2], v1, v2, v3);

            int idatom = state.raw.idfirst[itype - 1] - 1;
            for (int iatom = 1; iatom <= state.raw.natoms[itype - 1]; ++iatom) {
                icart += 1;
                idatom += 1;

                const std::array<double, 3> xref = {
                    coor_value(state, idatom, 1),
                    coor_value(state, idatom, 2),
                    coor_value(state, idatom, 3),
                };
                double out[3] = {0.0, 0.0, 0.0};
                compcart_cpp(out, xcm, xref, v1, v2, v3);
                for (int axis = 1; axis <= 3; ++axis) {
                    xcart_ref(state, icart, axis) = out[axis - 1];
                }

                const double fplus = comprest_atom_cpp(state, icart);
                *f += fplus;
                packmol_frest = std::max(packmol_frest, fplus);
                packmol_constraint_penalty_sum += fplus;
                if (fplus > 0.0) {
                    packmol_constraint_penalty_count += 1;
                }
                if (state.meta.move && state.raw.frest_atom != nullptr) {
                    state.raw.frest_atom[icart - 1] += fplus;
                }

                if (!state.meta.init1) {
                    int cell[3] = {0, 0, 0};
                    setcell_cpp(state, out, cell);
                    const int latomfirst_idx = latomfirst_index(state, cell[0], cell[1], cell[2]);
                    const bool empty_cell = state.raw.latomfirst[latomfirst_idx] == 0;
                    state.raw.latomnext[icart - 1] = state.raw.latomfirst[latomfirst_idx];
                    state.raw.latomfirst[latomfirst_idx] = icart;
                    if (empty_cell) {
                        const int icell = index_cell_cpp(state, cell);
                        state.raw.lcellnext[icell - 1] = lcellfirst;
                        lcellfirst = icell;
                    }
                    state.raw.ibtype[icart - 1] = itype;
                    state.raw.ibmol[icart - 1] = imol;
                }
            }

            ilugan += 3;
            ilubar += 3;
        }
    }

    if (state.meta.init1) {
        return true;
    }

    int icell = lcellfirst;
    while (icell > 0) {
        int cell[3] = {0, 0, 0};
        icell_to_cell_cpp(state, icell, cell);
        const int i = cell[0];
        const int j = cell[1];
        const int k = cell[2];

        int icart_cell = state.raw.latomfirst[latomfirst_index(state, i, j, k)];
        while (icart_cell > 0) {
            if (state.comptype_flags[state.raw.ibtype[icart_cell - 1] - 1] != 0) {
                *f += fparc_cpp(state, icart_cell, state.raw.latomnext[icart_cell - 1]);
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), j, k
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state, i, cell_ind_cpp(j + 1, state.meta.ncells[1]), k
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state, i, j, cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        cell_ind_cpp(j - 1, state.meta.ncells[1]),
                        k
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        j,
                        cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        i,
                        cell_ind_cpp(j + 1, state.meta.ncells[1]),
                        cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        i,
                        cell_ind_cpp(j + 1, state.meta.ncells[1]),
                        cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        cell_ind_cpp(j + 1, state.meta.ncells[1]),
                        k
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        j,
                        cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        cell_ind_cpp(j - 1, state.meta.ncells[1]),
                        cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        cell_ind_cpp(j - 1, state.meta.ncells[1]),
                        cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        cell_ind_cpp(j + 1, state.meta.ncells[1]),
                        cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]
                );
                *f += fparc_cpp(
                    state,
                    icart_cell,
                    state.raw.latomfirst[latomfirst_index(
                        state,
                        cell_ind_cpp(i + 1, state.meta.ncells[0]),
                        cell_ind_cpp(j + 1, state.meta.ncells[1]),
                        cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]
                );
            }
            icart_cell = state.raw.latomnext[icart_cell - 1];
        }
        icell = state.raw.lcellnext[icell - 1];
    }

    packmol_eval_set_lcellfirst_c(&lcellfirst);

    return true;
}

} // namespace

bool try_eval_objective_full_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
) {
    if (!use_cpp_eval_kernel()) {
        return false;
    }
    return eval_objective_full_cpp_kernel(n, x, m, lambda, rho, f, inform);
}

bool eval_objective_full_cpp_direct(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* f,
    int* inform
) {
    return eval_objective_full_cpp_kernel(n, x, m, lambda, rho, f, inform);
}

namespace {

double& gxcar_ref(const PackmolEvalStateSnapshot& state, const int icart, const int axis) {
    return state.raw.gxcar[(icart - 1) + (axis - 1) * state.meta.ntotat];
}

void gwalls_cpp(const PackmolEvalStateSnapshot& state, const int icart, const int irest) {
    const int type = state.raw.ityperest[irest - 1];
    const double scale = state.meta.scale;
    const double scale2 = state.meta.scale2;
    const double x = xcart_value(state, icart, 1);
    const double y = xcart_value(state, icart, 2);
    const double z = xcart_value(state, icart, 3);

    if (type == 2 || type == 3) {
        const double clength = (type == 2) ? restpars_value(state, irest, 4) : 0.0;
        const double xmin = restpars_value(state, irest, 1);
        const double ymin = restpars_value(state, irest, 2);
        const double zmin = restpars_value(state, irest, 3);
        const double xmax = (type == 2) ? (xmin + clength) : restpars_value(state, irest, 4);
        const double ymax = (type == 2) ? (ymin + clength) : restpars_value(state, irest, 5);
        const double zmax = (type == 2) ? (zmin + clength) : restpars_value(state, irest, 6);

        double a1 = x - xmin;
        double a2 = y - ymin;
        double a3 = z - zmin;
        if (a1 < 0.0) {
            gxcar_ref(state, icart, 1) += scale * 2.0 * a1;
        }
        if (a2 < 0.0) {
            gxcar_ref(state, icart, 2) += scale * 2.0 * a2;
        }
        if (a3 < 0.0) {
            gxcar_ref(state, icart, 3) += scale * 2.0 * a3;
        }
        a1 = x - xmax;
        a2 = y - ymax;
        a3 = z - zmax;
        if (a1 > 0.0) {
            gxcar_ref(state, icart, 1) += scale * 2.0 * a1;
        }
        if (a2 > 0.0) {
            gxcar_ref(state, icart, 2) += scale * 2.0 * a2;
        }
        if (a3 > 0.0) {
            gxcar_ref(state, icart, 3) += scale * 2.0 * a3;
        }
        return;
    }

    if (type == 4 || type == 8) {
        const double d =
            std::pow(x - restpars_value(state, irest, 1), 2) +
            std::pow(y - restpars_value(state, irest, 2), 2) +
            std::pow(z - restpars_value(state, irest, 3), 2) -
            std::pow(restpars_value(state, irest, 4), 2);
        if ((type == 4 && d > 0.0) || (type == 8 && d < 0.0)) {
            gxcar_ref(state, icart, 1) += 4.0 * scale2 * (x - restpars_value(state, irest, 1)) * d;
            gxcar_ref(state, icart, 2) += 4.0 * scale2 * (y - restpars_value(state, irest, 2)) * d;
            gxcar_ref(state, icart, 3) += 4.0 * scale2 * (z - restpars_value(state, irest, 3)) * d;
        }
        return;
    }

    if (type == 5 || type == 9) {
        const double a1 = x - restpars_value(state, irest, 1);
        const double b1 = y - restpars_value(state, irest, 2);
        const double c1 = z - restpars_value(state, irest, 3);
        const double a2 = std::pow(restpars_value(state, irest, 4), 2);
        const double b2 = std::pow(restpars_value(state, irest, 5), 2);
        const double c2 = std::pow(restpars_value(state, irest, 6), 2);
        double d = a1 * a1 / a2 + b1 * b1 / b2 + c1 * c1 / c2 -
                   std::pow(restpars_value(state, irest, 7), 2);
        if ((type == 5 && d > 0.0) || (type == 9 && d < 0.0)) {
            if (type == 9) {
                d = scale2 * d;
                gxcar_ref(state, icart, 1) += 4.0 * d * a1 / a2;
                gxcar_ref(state, icart, 2) += 4.0 * d * b1 / b2;
                gxcar_ref(state, icart, 3) += 4.0 * d * c1 / c2;
            } else {
                gxcar_ref(state, icart, 1) += scale2 * 4.0 * d * a1 / a2;
                gxcar_ref(state, icart, 2) += scale2 * 4.0 * d * b1 / b2;
                gxcar_ref(state, icart, 3) += scale2 * 4.0 * d * c1 / c2;
            }
        }
        return;
    }

    if (type == 6 || type == 7) {
        const double xmin = restpars_value(state, irest, 1);
        const double ymin = restpars_value(state, irest, 2);
        const double zmin = restpars_value(state, irest, 3);
        const double xmax = (type == 6) ? (xmin + restpars_value(state, irest, 4)) : restpars_value(state, irest, 4);
        const double ymax = (type == 6) ? (ymin + restpars_value(state, irest, 4)) : restpars_value(state, irest, 5);
        const double zmax = (type == 6) ? (zmin + restpars_value(state, irest, 4)) : restpars_value(state, irest, 6);
        double a1 = 0.0;
        double a2 = 0.0;
        double a3 = 0.0;
        double a4 = 0.0;
        double a5 = 0.0;
        double a6 = 0.0;
        if ((x > xmin && x < xmax) && (y > ymin && y < ymax) && (z > zmin && z < zmax)) {
            const double xmed = (xmax - xmin) / 2.0;
            const double ymed = (ymax - ymin) / 2.0;
            const double zmed = (zmax - zmin) / 2.0;
            if (x <= xmed) {
                a1 = 1.0;
            } else {
                a4 = -1.0;
            }
            if (y <= ymed) {
                a2 = 1.0;
            } else {
                a5 = -1.0;
            }
            if (z <= zmed) {
                a3 = 1.0;
            } else {
                a6 = -1.0;
            }
        }
        gxcar_ref(state, icart, 1) += scale * (a1 + a4);
        gxcar_ref(state, icart, 2) += scale * (a2 + a5);
        gxcar_ref(state, icart, 3) += scale * (a3 + a6);
        return;
    }

    if (type == 10 || type == 11) {
        double d = restpars_value(state, irest, 1) * x +
                   restpars_value(state, irest, 2) * y +
                   restpars_value(state, irest, 3) * z -
                   restpars_value(state, irest, 4);
        if ((type == 10 && d < 0.0) || (type == 11 && d > 0.0)) {
            d = scale * d;
            gxcar_ref(state, icart, 1) += 2.0 * restpars_value(state, irest, 1) * d;
            gxcar_ref(state, icart, 2) += 2.0 * restpars_value(state, irest, 2) * d;
            gxcar_ref(state, icart, 3) += 2.0 * restpars_value(state, irest, 3) * d;
        }
        return;
    }

    if (type == 12 || type == 13) {
        const double a1 = x - restpars_value(state, irest, 1);
        const double a2 = y - restpars_value(state, irest, 2);
        const double a3 = z - restpars_value(state, irest, 3);
        const double vnorm = std::sqrt(
            std::pow(restpars_value(state, irest, 4), 2) +
            std::pow(restpars_value(state, irest, 5), 2) +
            std::pow(restpars_value(state, irest, 6), 2)
        );
        const double vv1 = restpars_value(state, irest, 4) / vnorm;
        const double vv2 = restpars_value(state, irest, 5) / vnorm;
        const double vv3 = restpars_value(state, irest, 6) / vnorm;
        const double b1 = vv1 * a1;
        const double b2 = vv2 * a2;
        const double b3 = vv3 * a3;
        const double w = b1 + b2 + b3;
        const double d = std::pow(a1 - vv1 * w, 2) + std::pow(a2 - vv2 * w, 2) +
                         std::pow(a3 - vv3 * w, 2);

        if (type == 12) {
            const double rg1 = scale2 * (
                -2.0 * std::max(-w, 0.0) * vv1 +
                2.0 * std::max(w - restpars_value(state, irest, 9), 0.0) * vv1 +
                2.0 * std::max(d - std::pow(restpars_value(state, irest, 7), 2), 0.0) *
                    (2.0 * (a1 - vv1 * w) * (1.0 - vv1 * vv1) +
                     2.0 * (a2 - vv2 * w) * (-vv2 * vv1) +
                     2.0 * (a3 - vv3 * w) * (-vv3 * vv1))
            );
            const double rg2 = scale2 * (
                -2.0 * std::max(-w, 0.0) * vv2 +
                2.0 * std::max(w - restpars_value(state, irest, 9), 0.0) * vv2 +
                2.0 * std::max(d - std::pow(restpars_value(state, irest, 7), 2), 0.0) *
                    (2.0 * (a1 - vv1 * w) * (-vv1 * vv2) +
                     2.0 * (a2 - vv2 * w) * (1.0 - vv2 * vv2) +
                     2.0 * (a3 - vv3 * w) * (-vv3 * vv2))
            );
            const double rg3 = scale2 * (
                -2.0 * std::max(-w, 0.0) * vv3 +
                2.0 * std::max(w - restpars_value(state, irest, 9), 0.0) * vv3 +
                2.0 * std::max(d - std::pow(restpars_value(state, irest, 7), 2), 0.0) *
                    (2.0 * (a1 - vv1 * w) * (-vv1 * vv3) +
                     2.0 * (a2 - vv2 * w) * (-vv2 * vv3) +
                     2.0 * (a3 - vv3 * w) * (1.0 - vv3 * vv3))
            );
            gxcar_ref(state, icart, 1) += rg1;
            gxcar_ref(state, icart, 2) += rg2;
            gxcar_ref(state, icart, 3) += rg3;
        } else {
            const double frab = std::pow(std::min(-w, 0.0), 2) *
                                std::pow(std::min(w - restpars_value(state, irest, 9), 0.0), 2);
            const double frac = std::pow(std::min(-w, 0.0), 2) *
                                std::pow(std::min(d - std::pow(restpars_value(state, irest, 7), 2), 0.0), 2);
            const double frbc = std::pow(std::min(w - restpars_value(state, irest, 9), 0.0), 2) *
                                std::pow(std::min(d - std::pow(restpars_value(state, irest, 7), 2), 0.0), 2);
            const double dfra1 = -2.0 * std::min(-w, 0.0) * vv1;
            const double dfrb1 = 2.0 * std::min(w - restpars_value(state, irest, 9), 0.0) * vv1;
            const double dfrc1 = 2.0 * std::min(d - std::pow(restpars_value(state, irest, 7), 2), 0.0) *
                                 (2.0 * (a1 - vv1 * w) * (1.0 - vv1 * vv1) +
                                  2.0 * (a2 - vv2 * w) * (-vv2 * vv1) +
                                  2.0 * (a3 - vv3 * w) * (-vv3 * vv1));
            const double dfra2 = -2.0 * std::min(-w, 0.0) * vv2;
            const double dfrb2 = 2.0 * std::min(w - restpars_value(state, irest, 9), 0.0) * vv2;
            const double dfrc2 = 2.0 * std::min(d - std::pow(restpars_value(state, irest, 7), 2), 0.0) *
                                 (2.0 * (a1 - vv1 * w) * (-vv1 * vv2) +
                                  2.0 * (a2 - vv2 * w) * (1.0 - vv2 * vv2) +
                                  2.0 * (a3 - vv3 * w) * (-vv3 * vv2));
            const double dfra3 = -2.0 * std::min(-w, 0.0) * vv3;
            const double dfrb3 = 2.0 * std::min(w - restpars_value(state, irest, 9), 0.0) * vv3;
            const double dfrc3 = 2.0 * std::min(d - std::pow(restpars_value(state, irest, 7), 2), 0.0) *
                                 (2.0 * (a1 - vv1 * w) * (-vv1 * vv3) +
                                  2.0 * (a2 - vv2 * w) * (-vv2 * vv3) +
                                  2.0 * (a3 - vv3 * w) * (1.0 - vv3 * vv3));
            gxcar_ref(state, icart, 1) += scale2 * (dfra1 * frbc + dfrb1 * frac + dfrc1 * frab);
            gxcar_ref(state, icart, 2) += scale2 * (dfra2 * frbc + dfrb2 * frac + dfrc2 * frab);
            gxcar_ref(state, icart, 3) += scale2 * (dfra3 * frbc + dfrb3 * frac + dfrc3 * frab);
        }
        return;
    }

    if (type == 14 || type == 15) {
        const double a1 = -std::pow(x - restpars_value(state, irest, 1), 2) /
                          (2.0 * std::pow(restpars_value(state, irest, 3), 2));
        const double a2 = -std::pow(y - restpars_value(state, irest, 2), 2) /
                          (2.0 * std::pow(restpars_value(state, irest, 4), 2));
        double d = 0.0;
        if (a1 + a2 <= -50.0) {
            d = -(z - restpars_value(state, irest, 5));
        } else {
            d = restpars_value(state, irest, 6) * std::exp(a1 + a2) -
                (z - restpars_value(state, irest, 5));
        }
        if ((type == 14 && d > 0.0) || (type == 15 && d < 0.0)) {
            d = scale * d;
            gxcar_ref(state, icart, 1) -= 2.0 * d * (x - restpars_value(state, irest, 1)) *
                                          (d + (z - restpars_value(state, irest, 5))) /
                                          std::pow(restpars_value(state, irest, 3), 2);
            gxcar_ref(state, icart, 2) -= 2.0 * d * (y - restpars_value(state, irest, 2)) *
                                          (d + (z - restpars_value(state, irest, 5))) /
                                          std::pow(restpars_value(state, irest, 4), 2);
            gxcar_ref(state, icart, 3) -= 2.0 * d;
        }
    }
}

void gparc_cpp_gradient(const PackmolEvalStateSnapshot& state, const int icart, int firstjcart) {
    int jcart = firstjcart;
    while (jcart > 0) {
        if (state.comptype_flags[state.raw.ibtype[jcart - 1] - 1] == 0) {
            jcart = state.raw.latomnext[jcart - 1];
            continue;
        }
        if (state.raw.ibmol[icart - 1] == state.raw.ibmol[jcart - 1] &&
            state.raw.ibtype[icart - 1] == state.raw.ibtype[jcart - 1]) {
            jcart = state.raw.latomnext[jcart - 1];
            continue;
        }
        if (state.fixedatom_flags[icart - 1] != 0 && state.fixedatom_flags[jcart - 1] != 0) {
            jcart = state.raw.latomnext[jcart - 1];
            continue;
        }

        const double dx = delta_vector_cpp(
            xcart_value(state, icart, 1), xcart_value(state, jcart, 1), state.meta.pbc_length[0]
        );
        const double dy = delta_vector_cpp(
            xcart_value(state, icart, 2), xcart_value(state, jcart, 2), state.meta.pbc_length[1]
        );
        const double dz = delta_vector_cpp(
            xcart_value(state, icart, 3), xcart_value(state, jcart, 3), state.meta.pbc_length[2]
        );
        const double datom = dx * dx + dy * dy + dz * dz;
        const double tol = std::pow(state.raw.radius[icart - 1] + state.raw.radius[jcart - 1], 2);
        if (datom < tol) {
            double dtemp = state.raw.fscale[icart - 1] * state.raw.fscale[jcart - 1] * 4.0 *
                           (datom - tol);
            double xdiff = dtemp * dx;
            gxcar_ref(state, icart, 1) += xdiff;
            gxcar_ref(state, jcart, 1) -= xdiff;
            xdiff = dtemp * dy;
            gxcar_ref(state, icart, 2) += xdiff;
            gxcar_ref(state, jcart, 2) -= xdiff;
            xdiff = dtemp * dz;
            gxcar_ref(state, icart, 3) += xdiff;
            gxcar_ref(state, jcart, 3) -= xdiff;

            if (state.use_short_radius_flags[icart - 1] != 0 ||
                state.use_short_radius_flags[jcart - 1] != 0) {
                const double short_tol = std::pow(
                    state.raw.short_radius[icart - 1] + state.raw.short_radius[jcart - 1], 2
                );
                if (datom < short_tol) {
                    double short_tol_scale = std::sqrt(
                        state.raw.short_radius_scale[icart - 1] *
                        state.raw.short_radius_scale[jcart - 1]
                    );
                    short_tol_scale *= (tol * tol) / (short_tol * short_tol);
                    dtemp = state.raw.fscale[icart - 1] * state.raw.fscale[jcart - 1] * 4.0 *
                            short_tol_scale * (datom - short_tol);
                    xdiff = dtemp * dx;
                    gxcar_ref(state, icart, 1) += xdiff;
                    gxcar_ref(state, jcart, 1) -= xdiff;
                    xdiff = dtemp * dy;
                    gxcar_ref(state, icart, 2) += xdiff;
                    gxcar_ref(state, jcart, 2) -= xdiff;
                    xdiff = dtemp * dz;
                    gxcar_ref(state, icart, 3) += xdiff;
                    gxcar_ref(state, jcart, 3) -= xdiff;
                }
            }
        }
        jcart = state.raw.latomnext[jcart - 1];
    }
}

bool eval_gradient_full_cpp_kernel(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    int* inform
) {
    (void)m;
    (void)lambda;
    (void)rho;

    if (n == nullptr || x == nullptr || g == nullptr || inform == nullptr) {
        return false;
    }

    PackmolEvalStateSnapshot state;
    if (!load_packmol_eval_state_cpp(&state)) {
        return false;
    }
    if (state.raw.gxcar == nullptr || *n != state.meta.ntotmol * 6) {
        return false;
    }

    *inform = 0;
    for (int i = 0; i < state.meta.ntotat * 3; ++i) {
        state.raw.gxcar[i] = 0.0;
    }

    int lcellfirst = state.meta.lcellfirst;
    if (!state.meta.init1) {
        resetcells_cpp(state, &lcellfirst);
    }

    int ilubar = 0;
    int ilugan = state.meta.ntotmol * 3;
    int icart = 0;
    for (int itype = 1; itype <= state.meta.ntype; ++itype) {
        if (state.comptype_flags[itype - 1] == 0) {
            icart += state.raw.nmols[itype - 1] * state.raw.natoms[itype - 1];
            continue;
        }

        for (int imol = 1; imol <= state.raw.nmols[itype - 1]; ++imol) {
            const std::array<double, 3> xcm = {x[ilubar], x[ilubar + 1], x[ilubar + 2]};
            std::array<double, 3> v1{};
            std::array<double, 3> v2{};
            std::array<double, 3> v3{};
            eulerrmat_cpp(x[ilugan], x[ilugan + 1], x[ilugan + 2], v1, v2, v3);
            int idatom = state.raw.idfirst[itype - 1] - 1;

            for (int iatom = 1; iatom <= state.raw.natoms[itype - 1]; ++iatom) {
                icart += 1;
                idatom += 1;
                const std::array<double, 3> xref = {
                    coor_value(state, idatom, 1),
                    coor_value(state, idatom, 2),
                    coor_value(state, idatom, 3),
                };
                double out[3] = {0.0, 0.0, 0.0};
                compcart_cpp(out, xcm, xref, v1, v2, v3);
                for (int axis = 1; axis <= 3; ++axis) {
                    xcart_ref(state, icart, axis) = out[axis - 1];
                }

                for (int iratcount = 1; iratcount <= state.raw.nratom[icart - 1]; ++iratcount) {
                    gwalls_cpp(state, icart, iratom_value(state, icart, iratcount));
                }

                if (!state.meta.init1) {
                    int cell[3] = {0, 0, 0};
                    setcell_cpp(state, out, cell);
                    const int latomfirst_idx = latomfirst_index(state, cell[0], cell[1], cell[2]);
                    const bool empty_cell = state.raw.latomfirst[latomfirst_idx] == 0;
                    state.raw.latomnext[icart - 1] = state.raw.latomfirst[latomfirst_idx];
                    state.raw.latomfirst[latomfirst_idx] = icart;
                    if (empty_cell) {
                        const int icell = index_cell_cpp(state, cell);
                        state.raw.lcellnext[icell - 1] = lcellfirst;
                        lcellfirst = icell;
                    }
                    state.raw.ibtype[icart - 1] = itype;
                    state.raw.ibmol[icart - 1] = imol;
                }
            }
            ilugan += 3;
            ilubar += 3;
        }
    }

    if (!state.meta.init1) {
        int icell = lcellfirst;
        while (icell > 0) {
            int cell[3] = {0, 0, 0};
            icell_to_cell_cpp(state, icell, cell);
            const int i = cell[0];
            const int j = cell[1];
            const int k = cell[2];

            int icart_cell = state.raw.latomfirst[latomfirst_index(state, i, j, k)];
            while (icart_cell > 0) {
                if (state.comptype_flags[state.raw.ibtype[icart_cell - 1] - 1] != 0) {
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomnext[icart_cell - 1]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), j, k
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, i, cell_ind_cpp(j + 1, state.meta.ncells[1]), k
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, i, j, cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, i, cell_ind_cpp(j + 1, state.meta.ncells[1]), cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, i, cell_ind_cpp(j + 1, state.meta.ncells[1]), cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), cell_ind_cpp(j + 1, state.meta.ncells[1]), k
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), j, cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), cell_ind_cpp(j - 1, state.meta.ncells[1]), k
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), j, cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), cell_ind_cpp(j + 1, state.meta.ncells[1]), cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), cell_ind_cpp(j + 1, state.meta.ncells[1]), cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), cell_ind_cpp(j - 1, state.meta.ncells[1]), cell_ind_cpp(k + 1, state.meta.ncells[2])
                    )]);
                    gparc_cpp_gradient(state, icart_cell, state.raw.latomfirst[latomfirst_index(
                        state, cell_ind_cpp(i + 1, state.meta.ncells[0]), cell_ind_cpp(j - 1, state.meta.ncells[1]), cell_ind_cpp(k - 1, state.meta.ncells[2])
                    )]);
                }
                icart_cell = state.raw.latomnext[icart_cell - 1];
            }
            icell = state.raw.lcellnext[icell - 1];
        }
    }

    std::fill(g, g + *n, 0.0);
    int k1 = 0;
    int k2 = state.meta.ntotmol * 3;
    icart = 0;
    for (int itype = 1; itype <= state.meta.ntype; ++itype) {
        if (state.comptype_flags[itype - 1] == 0) {
            icart += state.raw.nmols[itype - 1] * state.raw.natoms[itype - 1];
            continue;
        }

        for (int imol = 1; imol <= state.raw.nmols[itype - 1]; ++imol) {
            const double beta = x[k2];
            const double gama = x[k2 + 1];
            const double teta = x[k2 + 2];
            const double cb = std::cos(beta);
            const double sb = std::sin(beta);
            const double cg = std::cos(gama);
            const double sg = std::sin(gama);
            const double ct = std::cos(teta);
            const double st = std::sin(teta);

            double dv1beta[3] = {-cb * sg * ct - sb * cg, -cb * cg * ct + sb * sg, cb * st};
            double dv2beta[3] = {-sb * sg * ct + cb * cg, -sb * cg * ct - cb * sg, sb * st};
            double dv3beta[3] = {0.0, 0.0, 0.0};

            double dv1gama[3] = {-sb * cg * ct - cb * sg, sb * sg * ct - cb * cg, 0.0};
            double dv2gama[3] = {cb * cg * ct - sb * sg, -sg * cb * ct - cg * sb, 0.0};
            double dv3gama[3] = {cg * st, -sg * st, 0.0};

            double dv1teta[3] = {sb * sg * st, sb * cg * st, sb * ct};
            double dv2teta[3] = {-cb * sg * st, -cb * cg * st, -cb * ct};
            double dv3teta[3] = {sg * ct, cg * ct, -st};

            int idatom = state.raw.idfirst[itype - 1] - 1;
            for (int iatom = 1; iatom <= state.raw.natoms[itype - 1]; ++iatom) {
                icart += 1;
                idatom += 1;

                for (int k = 0; k < 3; ++k) {
                    g[k1 + k] += gxcar_ref(state, icart, k + 1);
                }

                for (int k = 0; k < 3; ++k) {
                    g[k2] += (coor_value(state, idatom, 1) * dv1beta[k] +
                              coor_value(state, idatom, 2) * dv2beta[k] +
                              coor_value(state, idatom, 3) * dv3beta[k]) *
                             gxcar_ref(state, icart, k + 1);
                    g[k2 + 1] += (coor_value(state, idatom, 1) * dv1gama[k] +
                                  coor_value(state, idatom, 2) * dv2gama[k] +
                                  coor_value(state, idatom, 3) * dv3gama[k]) *
                                 gxcar_ref(state, icart, k + 1);
                    g[k2 + 2] += (coor_value(state, idatom, 1) * dv1teta[k] +
                                  coor_value(state, idatom, 2) * dv2teta[k] +
                                  coor_value(state, idatom, 3) * dv3teta[k]) *
                                 gxcar_ref(state, icart, k + 1);
                }
            }
            k2 += 3;
            k1 += 3;
        }
    }

    packmol_eval_set_lcellfirst_c(&lcellfirst);
    return true;
}

} // namespace

bool try_eval_gradient_full_cpp(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    int* inform
) {
    if (!use_cpp_gradient_kernel()) {
        return false;
    }
    return eval_gradient_full_cpp_kernel(n, x, m, lambda, rho, g, inform);
}

bool eval_gradient_full_cpp_direct(
    const int* n,
    double* x,
    const int* m,
    const double* lambda,
    const double* rho,
    double* g,
    int* inform
) {
    return eval_gradient_full_cpp_kernel(n, x, m, lambda, rho, g, inform);
}
