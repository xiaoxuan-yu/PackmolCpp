#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import pathlib
import shutil
import subprocess
import tempfile


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
    nested_tn_cpp_tail: bool,
) -> tuple[dict[str, float | int], str]:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    env["PACKMOL_GENCAN_NUMERIC_CPP"] = "0"
    if nested_tn_cpp_tail:
        env["PACKMOL_GENCAN_DEBUG"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH"] = "1"
        env["PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY"] = "1"
        env["PACKMOL_GENCAN_SPG_POST_CPP_HANDOFF"] = "1"
        env["PACKMOL_GENCAN_HANDOFF_CPP_TAIL"] = "1"
        env["PACKMOL_GENCAN_HANDOFF_CPP_TAIL_ALLOW_TN_POST"] = "1"
        env["PACKMOL_GENCAN_SPG_POST_CPP_REPLAY_MAX_DEPTH"] = "2"

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
        raise RuntimeError(output)
    parsed = _parse_probe_output(output)
    required = ("inform", "f", "gpsupn", "xsum", "gnorm2", "iter", "fcnt", "gcnt", "cgcnt")
    missing = [k for k in required if k not in parsed]
    if missing:
        raise RuntimeError(f"missing keys: {missing}")
    return parsed, output


def _print_delta_row(key: str, baseline: float | int, candidate: float | int) -> None:
    if isinstance(baseline, int):
        delta = int(candidate) - int(baseline)
        status = "OK" if delta == 0 else "DIFF"
        print(f"{key:>6}: fortran={int(baseline):>6} nested_cpp={int(candidate):>6} delta={delta:>6} {status}")
        return

    baseline_f = float(baseline)
    candidate_f = float(candidate)
    delta = candidate_f - baseline_f
    status = "OK" if math.isclose(baseline_f, candidate_f, rel_tol=1e-10, abs_tol=1e-12) else "DIFF"
    print(
        f"{key:>6}: fortran={baseline_f:.16e} nested_cpp={candidate_f:.16e} "
        f"delta={delta:.16e} {status}"
    )


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    probe = repo / "build" / "tests" / "test_gencan_spg_post_ab_probe"
    fixture = repo / "tests" / "fixtures" / "tiny_packing.inp"
    structure = fixture.parent / "minimal_structure.pdb"

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        if structure.exists():
            shutil.copy2(structure, tmpdir / structure.name)

        baseline, _ = _run_probe(probe, staged_input, tmpdir, mode="fortran", nested_tn_cpp_tail=False)
        nested, nested_output = _run_probe(
            probe,
            staged_input,
            tmpdir,
            mode="cpp",
            nested_tn_cpp_tail=True,
        )

    if "[gencan-cpp-spg-post-cpp-replay]" not in nested_output:
        raise RuntimeError("missing spg_post cpp replay marker")
    if "[gencan-cpp-fortran-tail] reason=spg_post_nonterminal" in nested_output:
        raise RuntimeError("spg_post_nonterminal unexpectedly used fortran-tail")

    print("=== SPG post nested-TN cpp-tail report ===")
    print("This probe is informational: DIFF rows show the remaining gap before nested tn_post closure can be promoted.")
    for key in ("inform", "iter", "fcnt", "gcnt", "cgcnt"):
        _print_delta_row(key, int(baseline[key]), int(nested[key]))
    for key in ("f", "gpsupn", "xsum", "gnorm2"):
        _print_delta_row(key, float(baseline[key]), float(nested[key]))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
