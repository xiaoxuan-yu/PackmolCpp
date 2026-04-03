#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import tempfile


PROBES = (
    "test_gencan_ab_probe",
    "test_gencan_entry_stop_ab_probe",
    "test_gencan_entry_stop_inform2_ab_probe",
    "test_gencan_entry_stop_inform3_ab_probe",
)


def _parse_probe_output(text: str) -> dict[str, float | int]:
    values: dict[str, float | int] = {}
    for line in text.splitlines():
        tokens = line.strip().split()
        if not tokens:
            continue
        if tokens[0] == "mode" and len(tokens) >= 18:
            values["mode"] = int(tokens[1])
            values["inform"] = int(tokens[3])
            values["f"] = float(tokens[5])
            values["gpsupn"] = float(tokens[7])
            values["xsum"] = float(tokens[9])
            values["iter"] = int(tokens[11])
            values["fcnt"] = int(tokens[13])
            values["gcnt"] = int(tokens[15])
            values["cgcnt"] = int(tokens[17])
        elif tokens[0] == "gnorm2" and len(tokens) >= 2:
            values["gnorm2"] = float(tokens[1])
    return values


def _run_probe(
    probe: pathlib.Path,
    staged_input: pathlib.Path,
    workdir: pathlib.Path,
    *,
    mode: str,
    handoff: bool,
    canonicalize: bool,
    cpp_replay: bool = False,
    max_depth: int = 1,
    unsafe: bool = True,
) -> tuple[dict[str, float | int], str]:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    env["PACKMOL_GENCAN_NUMERIC_CPP"] = "0"
    env["PACKMOL_GENCAN_DEBUG"] = "1"
    if handoff:
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE"] = "1" if unsafe else "0"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH"] = str(max_depth)
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CANONICALIZE"] = "1" if canonicalize else "0"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY"] = "1" if cpp_replay else "0"
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
        raise RuntimeError(
            f"probe failed: mode={mode} handoff={handoff} canonicalize={canonicalize}\n"
            f"{output}"
        )
    parsed = _parse_probe_output(output)
    required = ("inform", "f", "gpsupn", "xsum", "gnorm2", "iter", "fcnt", "gcnt", "cgcnt")
    missing = [k for k in required if k not in parsed]
    if missing:
        raise RuntimeError(f"missing keys: {missing}")
    return parsed, output


def _fmt_float(v: float) -> str:
    return f"{v:.9e}"


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    fixture = repo / "tests" / "fixtures" / "tiny_packing.inp"
    probe_dir = repo / "build" / "tests"

    print("=== TN handoff matrix (fortran baseline) ===")
    for probe_name in PROBES:
        probe = probe_dir / probe_name
        with tempfile.TemporaryDirectory() as tmp_name:
            tmpdir = pathlib.Path(tmp_name)
            staged_input = tmpdir / fixture.name
            shutil.copy2(fixture, staged_input)
            structure_src = fixture.parent / "minimal_structure.pdb"
            if structure_src.exists():
                shutil.copy2(structure_src, tmpdir / structure_src.name)

            baseline, _ = _run_probe(
                probe, staged_input, tmpdir, mode="fortran", handoff=False, canonicalize=True
            )
            canon, _ = _run_probe(
                probe, staged_input, tmpdir, mode="cpp", handoff=True, canonicalize=True
            )
            replay, _ = _run_probe(
                probe,
                staged_input,
                tmpdir,
                mode="cpp",
                handoff=True,
                canonicalize=True,
                cpp_replay=True,
            )
            nocanon, _ = _run_probe(
                probe, staged_input, tmpdir, mode="cpp", handoff=True, canonicalize=False
            )
            safe_no_legacy, _ = _run_probe(
                probe,
                staged_input,
                tmpdir,
                mode="cpp",
                handoff=True,
                canonicalize=True,
                unsafe=False,
            )
            unsafe_no_legacy, unsafe_no_legacy_output = _run_probe(
                probe,
                staged_input,
                tmpdir,
                mode="cpp",
                handoff=True,
                canonicalize=True,
                unsafe=True,
            )

        if "[gencan-cpp-fortran-tail]" in unsafe_no_legacy_output:
            raise RuntimeError(f"{probe_name}: unexpected stale fortran-tail marker")
        for key in ("inform", "iter", "fcnt", "gcnt", "cgcnt"):
            if int(safe_no_legacy[key]) != int(unsafe_no_legacy[key]):
                raise RuntimeError(
                    f"{probe_name}: safe/unsafe drift in {key}: "
                    f"safe={safe_no_legacy[key]} unsafe={unsafe_no_legacy[key]}"
                )
        for key in ("f", "gpsupn", "xsum", "gnorm2"):
            lhs = float(safe_no_legacy[key])
            rhs = float(unsafe_no_legacy[key])
            if abs(lhs - rhs) > 1.0e-12:
                raise RuntimeError(
                    f"{probe_name}: safe/unsafe drift in {key}: "
                    f"safe={lhs:.16e} unsafe={rhs:.16e}"
                )

        print(f"\n[{probe_name}]")
        print(
            "  inform:"
            f" baseline={int(baseline['inform'])}"
            f" canon={int(canon['inform'])}"
            f" replay={int(replay['inform'])}"
            f" nocanon={int(nocanon['inform'])}"
        )
        print(
            "  iter/fcnt/gcnt/cgcnt:"
            f" baseline={int(baseline['iter'])}/{int(baseline['fcnt'])}/{int(baseline['gcnt'])}/{int(baseline['cgcnt'])}"
            f" canon={int(canon['iter'])}/{int(canon['fcnt'])}/{int(canon['gcnt'])}/{int(canon['cgcnt'])}"
            f" replay={int(replay['iter'])}/{int(replay['fcnt'])}/{int(replay['gcnt'])}/{int(replay['cgcnt'])}"
            f" nocanon={int(nocanon['iter'])}/{int(nocanon['fcnt'])}/{int(nocanon['gcnt'])}/{int(nocanon['cgcnt'])}"
        )
        print(
            "  f:"
            f" baseline={_fmt_float(float(baseline['f']))}"
            f" canon={_fmt_float(float(canon['f']))}"
            f" replay={_fmt_float(float(replay['f']))}"
            f" nocanon={_fmt_float(float(nocanon['f']))}"
        )
        print(
            "  gpsupn:"
            f" baseline={_fmt_float(float(baseline['gpsupn']))}"
            f" canon={_fmt_float(float(canon['gpsupn']))}"
            f" replay={_fmt_float(float(replay['gpsupn']))}"
            f" nocanon={_fmt_float(float(nocanon['gpsupn']))}"
        )
        print(
            "  xsum:"
            f" baseline={_fmt_float(float(baseline['xsum']))}"
            f" canon={_fmt_float(float(canon['xsum']))}"
            f" replay={_fmt_float(float(replay['xsum']))}"
            f" nocanon={_fmt_float(float(nocanon['xsum']))}"
        )
        print(
            "  gnorm2:"
            f" baseline={_fmt_float(float(baseline['gnorm2']))}"
            f" canon={_fmt_float(float(canon['gnorm2']))}"
            f" replay={_fmt_float(float(replay['gnorm2']))}"
            f" nocanon={_fmt_float(float(nocanon['gnorm2']))}"
        )
        print(
            "  legacy_off (safe/unsafe):"
            f" inform={int(safe_legacy_off['inform'])}/{int(unsafe_legacy_off['inform'])}"
            f" iter={int(safe_legacy_off['iter'])}/{int(unsafe_legacy_off['iter'])}"
            f" f={_fmt_float(float(safe_legacy_off['f']))}/{_fmt_float(float(unsafe_legacy_off['f']))}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
