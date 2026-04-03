#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile


MODE_ID = {"fortran": 0, "cpp": 1, "ab": 2}


def _run_probe(
    probe: pathlib.Path,
    input_path: pathlib.Path,
    mode: str,
    workdir: pathlib.Path,
) -> str:
    env = os.environ.copy()
    env["PACKMOL_GENCAN_IMPL"] = mode
    env["PACKMOL_GENCAN_DEBUG"] = "1"
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
    mode_marker = f"mode {MODE_ID[mode]} "
    if mode_marker not in output:
        print(output, file=sys.stderr)
        raise RuntimeError(f"unexpected mode marker for mode={mode}")
    return output


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

        # Baseline still runs for parity workflow consistency.
        _run_probe(probe, staged_input, "fortran", tmpdir)
        cpp_output = _run_probe(probe, staged_input, "cpp", tmpdir)
        ab_output = _run_probe(probe, staged_input, "ab", tmpdir)

        marker = "[gencan-cpp-fallback]"
        if marker in cpp_output:
            raise RuntimeError("fallback marker observed in cpp mode")
        if marker in ab_output:
            raise RuntimeError("fallback marker observed in ab mode")
        missing_reason_marker = "[gencan-cpp-missing-tail-reason]"
        if missing_reason_marker in cpp_output:
            raise RuntimeError("missing tail reason marker observed in cpp mode")
        if missing_reason_marker in ab_output:
            raise RuntimeError("missing tail reason marker observed in ab mode")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
