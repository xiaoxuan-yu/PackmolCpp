#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile


MODE_ID = {
    "fortran": 0,
    "cpp": 1,
}

INT_RE = re.compile(r"^[+-]?\d+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare GENCAN probe behavior between the Fortran and C++ implementations."
    )
    parser.add_argument("probe", type=pathlib.Path)
    parser.add_argument(
        "--fixture",
        type=pathlib.Path,
        default=None,
        help="Optional input fixture staged into a temporary working directory before each run.",
    )
    parser.add_argument("--int", dest="int_keys", action="append", default=[])
    parser.add_argument("--float", dest="float_specs", action="append", default=[])
    return parser.parse_args()


def parse_probe_output(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in text.replace("\x00", "").splitlines():
        tokens = raw_line.strip().split()
        if not tokens or len(tokens) % 2 != 0:
            continue
        for index in range(0, len(tokens), 2):
            values[tokens[index]] = tokens[index + 1]
    return values


def parse_float_spec(spec: str) -> tuple[str, float, float]:
    parts = spec.split(":")
    if len(parts) == 1:
        return parts[0], 1.0e-12, 1.0e-10
    if len(parts) == 2:
        return parts[0], float(parts[1]), 1.0e-10
    if len(parts) == 3:
        return parts[0], float(parts[1]), float(parts[2])
    raise ValueError(f"invalid float comparison spec: {spec}")


def stage_fixture_dir(fixture: pathlib.Path, workdir: pathlib.Path) -> pathlib.Path:
    for entry in fixture.parent.iterdir():
        if entry.is_file():
            shutil.copy2(entry, workdir / entry.name)
    return workdir / fixture.name


def run_probe(
    probe: pathlib.Path,
    mode: str,
    fixture: pathlib.Path | None,
) -> tuple[dict[str, str], str]:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode

    with tempfile.TemporaryDirectory() as tmp_name:
        workdir = pathlib.Path(tmp_name)
        command = [str(probe)]
        if fixture is not None:
            staged_fixture = stage_fixture_dir(fixture, workdir)
            command.append(str(staged_fixture))

        completed = subprocess.run(
            command,
            cwd=workdir,
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        output = (completed.stdout + completed.stderr).replace("\x00", "")
        if completed.returncode != 0:
            raise RuntimeError(f"{probe.name} failed in mode={mode}\n{output}")

    parsed = parse_probe_output(output)
    expected_mode = MODE_ID[mode]
    raw_mode = parsed.get("mode")
    if raw_mode is None or not INT_RE.match(raw_mode) or int(raw_mode) != expected_mode:
        raise RuntimeError(
            f"{probe.name} returned unexpected mode marker in mode={mode}\n"
            f"expected {expected_mode}, got {raw_mode!r}\n{output}"
        )
    return parsed, output


def require_key(values: dict[str, str], key: str, mode: str, output: str) -> str:
    if key not in values:
        raise RuntimeError(f"missing key {key!r} in mode={mode}\n{output}")
    return values[key]


def compare_int_keys(
    probe: pathlib.Path,
    int_keys: list[str],
    baseline: dict[str, str],
    candidate: dict[str, str],
    baseline_output: str,
    candidate_output: str,
) -> None:
    for key in int_keys:
        lhs_raw = require_key(baseline, key, "fortran", baseline_output)
        rhs_raw = require_key(candidate, key, "cpp", candidate_output)
        try:
            lhs = int(lhs_raw)
            rhs = int(rhs_raw)
        except ValueError as exc:
            raise RuntimeError(f"{probe.name}: {key} is not an integer field") from exc
        if lhs != rhs:
            raise RuntimeError(f"{probe.name}: {key} mismatch: fortran={lhs} cpp={rhs}")


def compare_float_keys(
    probe: pathlib.Path,
    float_specs: list[str],
    baseline: dict[str, str],
    candidate: dict[str, str],
    baseline_output: str,
    candidate_output: str,
) -> None:
    for spec in float_specs:
        key, atol, rtol = parse_float_spec(spec)
        lhs_raw = require_key(baseline, key, "fortran", baseline_output)
        rhs_raw = require_key(candidate, key, "cpp", candidate_output)
        lhs = float(lhs_raw)
        rhs = float(rhs_raw)
        if not math.isfinite(lhs) or not math.isfinite(rhs):
            raise RuntimeError(f"{probe.name}: {key} must stay finite: fortran={lhs} cpp={rhs}")
        if not math.isclose(lhs, rhs, abs_tol=atol, rel_tol=rtol):
            raise RuntimeError(
                f"{probe.name}: {key} mismatch: fortran={lhs:.16e} cpp={rhs:.16e} "
                f"(atol={atol:.1e}, rtol={rtol:.1e})"
            )


def main() -> int:
    args = parse_args()
    probe = args.probe.resolve()
    fixture = args.fixture.resolve() if args.fixture is not None else None

    baseline, baseline_output = run_probe(probe, "fortran", fixture)
    candidate, candidate_output = run_probe(probe, "cpp", fixture)

    compare_int_keys(probe, args.int_keys, baseline, candidate, baseline_output, candidate_output)
    compare_float_keys(probe, args.float_specs, baseline, candidate, baseline_output, candidate_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
