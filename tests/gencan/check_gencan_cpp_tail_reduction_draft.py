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
    "test_gencan_entry_stop_gtype2_ab_probe",
    "test_gencan_entry_stop_inform2_ab_probe",
    "test_gencan_entry_stop_inform3_ab_probe",
    "test_gencan_runtime_driver_probe",
    "test_gencan_spg_post_ab_probe",
)


def _run_matrix(
    probe_dir: pathlib.Path,
    staged_input: pathlib.Path,
    tmpdir: pathlib.Path,
    env: dict[str, str],
) -> tuple[bool, bool]:
    saw_fortran_tail = False
    saw_blocked_tail = False
    for probe_name in PROBES:
        completed = subprocess.run(
            [str(probe_dir / probe_name), str(staged_input)],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        output = (completed.stdout + completed.stderr).replace("\x00", "")
        if completed.returncode != 0:
            raise RuntimeError(f"{probe_name} failed\n{output}")
        if "[gencan-cpp-fortran-tail]" in output:
            saw_fortran_tail = True
        if "[gencan-cpp-tail-blocked]" in output:
            saw_blocked_tail = True
    return saw_fortran_tail, saw_blocked_tail


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    probe_dir = repo / "build" / "tests"
    fixture = repo / "tests" / "fixtures" / "tiny_packing.inp"
    structure = fixture.parent / "minimal_structure.pdb"

    with tempfile.TemporaryDirectory() as tmp_name:
        tmpdir = pathlib.Path(tmp_name)
        staged_input = tmpdir / fixture.name
        shutil.copy2(fixture, staged_input)
        if structure.exists():
            shutil.copy2(structure, tmpdir / structure.name)

        env = os.environ.copy()
        env.update(
            {
                "PACKMOL_GENCAN_IMPL": "cpp",
                "PACKMOL_GENCAN_NUMERIC_CPP": "0",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF": "1",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_SAFE": "1",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE": "0",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH": "1",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CANONICALIZE": "1",
                "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY": "1",
                "PACKMOL_GENCAN_DEBUG": "1",
            }
        )

        default_fortran_tail, default_blocked_tail = _run_matrix(probe_dir, staged_input, tmpdir, env)
        if default_fortran_tail:
            raise RuntimeError("unexpected fortran-tail marker with default cpp tail reduction")
        if not default_blocked_tail:
            raise RuntimeError("default cpp tail reduction run did not trigger any blocked-tail marker")

        env_off = dict(env)
        env_off["PACKMOL_GENCAN_CPP_TAIL_REDUCTION"] = "0"
        off_fortran_tail, _ = _run_matrix(probe_dir, staged_input, tmpdir, env_off)
        if off_fortran_tail and default_fortran_tail:
            raise RuntimeError("unexpected fortran-tail marker in default cpp tail reduction run")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
