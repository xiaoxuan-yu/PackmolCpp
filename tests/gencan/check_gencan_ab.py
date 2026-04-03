#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile


MODE_ID = {
    "fortran": 0,
    "cpp": 1,
    "ab": 2,
}


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


def _run_probe(probe: pathlib.Path, input_path: pathlib.Path, mode: str, workdir: pathlib.Path) -> dict[str, float | int]:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    # Keep AB parity checks pinned to legacy-compatible numeric accumulation mode.
    env["PACKMOL_GENCAN_NUMERIC_CPP"] = "0"
    if mode != "fortran":
        env["PACKMOL_GENCAN_AB_COMPARE"] = "1"
    completed = subprocess.run(
        [str(probe), str(input_path)],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=env,
    )
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        print(f"probe failed in mode={mode}", file=sys.stderr)
        print(output, file=sys.stderr)
        raise RuntimeError("probe run failed")
    parsed = _parse_probe_output(output)
    expected_mode = MODE_ID[mode]
    if parsed.get("mode") != expected_mode:
        print(output, file=sys.stderr)
        raise RuntimeError(f"unexpected mode id for mode={mode}: {parsed.get('mode')} != {expected_mode}")
    required = ("inform", "f", "gpsupn", "xsum", "gnorm2", "iter", "fcnt", "gcnt", "cgcnt")
    missing = [key for key in required if key not in parsed]
    if missing:
        print(output, file=sys.stderr)
        raise RuntimeError(f"missing probe keys in mode={mode}: {missing}")
    return parsed


def _assert_close(name: str, lhs: float, rhs: float, atol: float = 1.0e-12, rtol: float = 1.0e-10) -> None:
    if not math.isclose(lhs, rhs, rel_tol=rtol, abs_tol=atol):
        raise RuntimeError(f"{name} mismatch: {lhs} vs {rhs}")


def main() -> int:
    probe = pathlib.Path(sys.argv[1])
    fixture = pathlib.Path(sys.argv[2])

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        structure_src = fixture.parent / "minimal_structure.pdb"
        if structure_src.exists():
            shutil.copy2(structure_src, tmpdir / structure_src.name)

        baseline = _run_probe(probe, staged_input, "fortran", tmpdir)
        candidate = _run_probe(probe, staged_input, "cpp", tmpdir)
        ab_mode = _run_probe(probe, staged_input, "ab", tmpdir)

        for key in ("inform", "iter", "fcnt", "gcnt", "cgcnt"):
            if int(baseline[key]) != int(candidate[key]):
                raise RuntimeError(f"{key} mismatch between fortran and cpp: {baseline[key]} vs {candidate[key]}")
            if int(baseline[key]) != int(ab_mode[key]):
                raise RuntimeError(f"{key} mismatch between fortran and ab: {baseline[key]} vs {ab_mode[key]}")

        for key in ("f", "gpsupn", "xsum", "gnorm2"):
            _assert_close(key, float(baseline[key]), float(candidate[key]))
            _assert_close(f"{key}_ab", float(baseline[key]), float(ab_mode[key]))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
