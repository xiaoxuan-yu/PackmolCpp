#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import shutil
import subprocess
import tempfile
import math


PROBES = (
    "test_gencan_ab_probe",
    "test_gencan_entry_stop_ab_probe",
    "test_gencan_entry_stop_gtype2_ab_probe",
    "test_gencan_entry_stop_inform2_ab_probe",
    "test_gencan_entry_stop_inform3_ab_probe",
    "test_gencan_runtime_driver_probe",
    "test_gencan_spg_post_ab_probe",
)

AB_COMPARE_PROBES = {
    "test_gencan_ab_probe",
    "test_gencan_entry_stop_ab_probe",
    "test_gencan_entry_stop_gtype2_ab_probe",
    "test_gencan_entry_stop_inform2_ab_probe",
    "test_gencan_entry_stop_inform3_ab_probe",
    "test_gencan_runtime_driver_probe",
    "test_gencan_spg_post_ab_probe",
}

STRICT_COUNTER_PROBES = {
    "test_gencan_ab_probe",
    "test_gencan_entry_stop_ab_probe",
    "test_gencan_entry_stop_gtype2_ab_probe",
    "test_gencan_entry_stop_inform2_ab_probe",
    "test_gencan_entry_stop_inform3_ab_probe",
    "test_gencan_runtime_driver_probe",
}

COUNTER_ABS_TOL: dict[str, dict[str, int]] = {
    "test_gencan_spg_post_ab_probe": {
        "iter": 4,
        "fcnt": 8,
        "gcnt": 6,
        "cgcnt": 3,
    },
}

FLOAT_ABS_TOL: dict[str, dict[str, float]] = {
    "test_gencan_spg_post_ab_probe": {
        "f": 5.0e-2,
        "gpsupn": 5.0e-2,
        "xsum": 5.0e-2,
        "gnorm2": 5.0e-2,
    },
}

SUMMARY_RE = re.compile(
    r"mode\s+(?P<mode>-?\d+)\s+"
    r"inform\s+(?P<inform>-?\d+)\s+"
    r"f\s+(?P<f>[-+0-9.Ee]+)\s+"
    r"gpsupn\s+(?P<gpsupn>[-+0-9.Ee]+)\s+"
    r"xsum\s+(?P<xsum>[-+0-9.Ee]+)\s+"
    r"iter\s+(?P<iter>-?\d+)\s+"
    r"fcnt\s+(?P<fcnt>-?\d+)\s+"
    r"gcnt\s+(?P<gcnt>-?\d+)\s+"
    r"cgcnt\s+(?P<cgcnt>-?\d+)"
)

GNORM2_RE = re.compile(r"gnorm2\s+(?P<gnorm2>[-+0-9.Ee]+)")


def _run_probe(
    probe: pathlib.Path,
    staged_input: pathlib.Path,
    workdir: pathlib.Path,
    *,
    impl_mode: str,
) -> str:
    env = os.environ.copy()
    env.update({
        "PACKMOL_GENCAN_IMPL": impl_mode,
        "PACKMOL_GENCAN_NUMERIC_CPP": "0",
        "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF": "1",
        "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_SAFE": "1",
        "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE": "0",
        "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH": "1",
        "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CANONICALIZE": "1",
        "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY": "1",
        "PACKMOL_GENCAN_DEBUG": "1",
    })
    completed = subprocess.run(
        [str(probe), str(staged_input)],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    output = (completed.stdout + completed.stderr).replace("\x00", "")
    if completed.returncode != 0:
        raise RuntimeError(f"{probe.name} failed\n{output}")
    return output


def _parse_summary(output: str) -> dict[str, float | int] | None:
    summary_match = SUMMARY_RE.search(output)
    if summary_match is None:
        return None

    out: dict[str, float | int] = {
        "mode": int(summary_match.group("mode")),
        "inform": int(summary_match.group("inform")),
        "f": float(summary_match.group("f")),
        "gpsupn": float(summary_match.group("gpsupn")),
        "xsum": float(summary_match.group("xsum")),
        "iter": int(summary_match.group("iter")),
        "fcnt": int(summary_match.group("fcnt")),
        "gcnt": int(summary_match.group("gcnt")),
        "cgcnt": int(summary_match.group("cgcnt")),
    }
    gnorm2_match = GNORM2_RE.search(output)
    if gnorm2_match is not None:
        out["gnorm2"] = float(gnorm2_match.group("gnorm2"))
    return out


def _assert_summary_sane(probe_name: str, output: str) -> dict[str, float | int] | None:
    summary = _parse_summary(output)
    if summary is None:
        return None

    mode = int(summary["mode"])
    inform = int(summary["inform"])
    iter_count = int(summary["iter"])
    fcnt = int(summary["fcnt"])
    gcnt = int(summary["gcnt"])
    cgcnt = int(summary["cgcnt"])

    if mode != 1:
        raise RuntimeError(f"{probe_name}: unexpected mode={mode}, expected cpp mode 1")
    if inform < 0 or inform > 8:
        raise RuntimeError(f"{probe_name}: unexpected inform={inform}, expected [0,8]")
    if iter_count < 0 or fcnt < 0 or gcnt < 0 or cgcnt < 0:
        raise RuntimeError(
            f"{probe_name}: negative counters iter={iter_count} fcnt={fcnt} gcnt={gcnt} cgcnt={cgcnt}"
        )

    for key in ("f", "gpsupn", "xsum"):
        value = float(summary[key])
        if not math.isfinite(value):
            raise RuntimeError(f"{probe_name}: non-finite {key}={value}")

    if "gnorm2" in summary:
        gnorm2 = float(summary["gnorm2"])
        if not math.isfinite(gnorm2):
            raise RuntimeError(f"{probe_name}: non-finite gnorm2={gnorm2}")
    return summary


def _format_summary(summary: dict[str, float | int]) -> str:
    keys = ("mode", "inform", "iter", "fcnt", "gcnt", "cgcnt", "f", "gpsupn", "xsum", "gnorm2")
    parts: list[str] = []
    for key in keys:
        if key not in summary:
            continue
        value = summary[key]
        if isinstance(value, int):
            parts.append(f"{key}={value}")
        else:
            parts.append(f"{key}={value:.16e}")
    return " ".join(parts)


def _raise_ab_drift(
    probe_name: str,
    reason: str,
    cpp_summary: dict[str, float | int],
    fortran_summary: dict[str, float | int],
) -> None:
    raise RuntimeError(
        f"{probe_name}: {reason}\n"
        f"  cpp: {_format_summary(cpp_summary)}\n"
        f"  fortran: {_format_summary(fortran_summary)}"
    )


def _assert_ab_drift(
    probe_name: str,
    cpp_summary: dict[str, float | int],
    fortran_summary: dict[str, float | int],
) -> None:

    if int(cpp_summary["mode"]) != 1:
        _raise_ab_drift(
            probe_name,
            f"cpp mode drift, got mode={cpp_summary['mode']}",
            cpp_summary,
            fortran_summary,
        )
    if int(fortran_summary["mode"]) != 0:
        _raise_ab_drift(
            probe_name,
            f"fortran mode drift, got mode={fortran_summary['mode']}",
            cpp_summary,
            fortran_summary,
        )

    if int(cpp_summary["inform"]) != int(fortran_summary["inform"]):
        _raise_ab_drift(
            probe_name,
            f"AB drift in inform: cpp={cpp_summary['inform']} fortran={fortran_summary['inform']}",
            cpp_summary,
            fortran_summary,
        )

    strict_counter_compare = probe_name in STRICT_COUNTER_PROBES
    if strict_counter_compare:
        for key in ("iter", "fcnt", "gcnt", "cgcnt"):
            if int(cpp_summary[key]) != int(fortran_summary[key]):
                _raise_ab_drift(
                    probe_name,
                    f"AB drift in {key}: cpp={cpp_summary[key]} fortran={fortran_summary[key]}",
                    cpp_summary,
                    fortran_summary,
                )
    else:
        counter_tol = COUNTER_ABS_TOL.get(probe_name, {})
        for key, tol in counter_tol.items():
            diff = abs(int(cpp_summary[key]) - int(fortran_summary[key]))
            if diff > tol:
                _raise_ab_drift(
                    probe_name,
                    f"AB drift in {key}: |cpp-fortran|={diff} exceeds tol={tol}",
                    cpp_summary,
                    fortran_summary,
                )
    for key in ("f", "gpsupn", "xsum", "gnorm2"):
        if key not in cpp_summary or key not in fortran_summary:
            continue
        cpp_v = float(cpp_summary[key])
        ftn_v = float(fortran_summary[key])
        if not math.isfinite(cpp_v) or not math.isfinite(ftn_v):
            _raise_ab_drift(
                probe_name,
                f"AB non-finite {key}: cpp={cpp_v:.16e} fortran={ftn_v:.16e}",
                cpp_summary,
                fortran_summary,
            )
        float_tol = FLOAT_ABS_TOL.get(probe_name, {}).get(key)
        if float_tol is not None and abs(cpp_v - ftn_v) > float_tol:
            _raise_ab_drift(
                probe_name,
                f"AB drift in {key}: |cpp-fortran|={abs(cpp_v-ftn_v):.16e} exceeds tol={float_tol:.16e}",
                cpp_summary,
                fortran_summary,
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-out", default="", help="Write AB compare summaries as JSON.")
    args = parser.parse_args()

    repo = pathlib.Path(__file__).resolve().parents[2]
    probe_dir = repo / "build" / "tests"
    fixture = repo / "tests" / "fixtures" / "tiny_packing.inp"
    report_path = pathlib.Path(args.report_out) if args.report_out else None
    report_records: list[dict[str, object]] = []

    ab_compare_enabled = os.environ.get(
        "PACKMOL_GENCAN_AB_COMPARE",
        "0",
    ) in {"1", "true", "TRUE", "on", "ON", "yes", "YES"}
    ab_report_enabled = os.environ.get(
        "PACKMOL_GENCAN_AB_COMPARE_REPORT",
        "0",
    ) in {"1", "true", "TRUE", "on", "ON", "yes", "YES"}

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        structure = fixture.parent / "minimal_structure.pdb"
        if structure.exists():
            shutil.copy2(structure, tmpdir / structure.name)

        for probe_name in PROBES:
            output = _run_probe(
                probe_dir / probe_name,
                staged_input,
                tmpdir,
                impl_mode="cpp",
            )
            if "[gencan-cpp-fallback]" in output:
                raise RuntimeError(f"{probe_name}: unexpected fallback marker")
            if "[gencan-cpp-fallback-blocked]" in output:
                raise RuntimeError(f"{probe_name}: unexpected blocked-fallback marker")
            if "[gencan-cpp-missing-tail-reason]" in output:
                raise RuntimeError(f"{probe_name}: missing tail reason marker")
            if "[gencan-cpp-fortran-tail]" in output:
                raise RuntimeError(f"{probe_name}: stale fortran-tail marker observed")
            if "[gencan-cpp-handoff-canonicalize-fortran]" in output:
                raise RuntimeError(f"{probe_name}: unexpected fortran canonicalize marker")
            cpp_summary = _assert_summary_sane(probe_name, output)
            if ab_compare_enabled and probe_name in AB_COMPARE_PROBES:
                fortran_output = _run_probe(
                    probe_dir / probe_name,
                    staged_input,
                    tmpdir,
                    impl_mode="fortran",
                )
                fortran_summary = _parse_summary(fortran_output)
                if cpp_summary is None or fortran_summary is None:
                    raise RuntimeError(f"{probe_name}: AB compare enabled but summary line missing")
                _assert_ab_drift(probe_name, cpp_summary, fortran_summary)
                report_records.append({
                    "probe": probe_name,
                    "cpp": cpp_summary,
                    "fortran": fortran_summary,
                })
                if ab_report_enabled:
                    print(
                        f"[ab-drift-ok] {probe_name} "
                        f"cpp=({_format_summary(cpp_summary)}) "
                        f"fortran=({_format_summary(fortran_summary)})"
                    )

    if report_path is not None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report_records, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
