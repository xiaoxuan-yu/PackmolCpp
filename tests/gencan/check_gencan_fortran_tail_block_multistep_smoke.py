#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
import tempfile


SUMMARY_RE = re.compile(r"mode\s+(?P<mode>-?\d+)\s+inform\s+(?P<inform>-?\d+)\s+")


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
                "PACKMOL_GENCAN_BLOCK_CPP_TAIL": "1",
                "PACKMOL_GENCAN_BLOCK_TAIL_STEPS": "4",
                "PACKMOL_GENCAN_DEBUG": "1",
            }
        )
        completed = subprocess.run(
            [str(probe), str(staged_input)],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        output = (completed.stdout + completed.stderr).replace("\x00", "")
        if "[gencan-cpp-tail-blocked]" not in output:
            raise RuntimeError("missing cpp-tail-blocked marker")
        if "[gencan-cpp-fortran-tail]" in output:
            raise RuntimeError("unexpected stale fortran-tail marker while blocked")

        summary = SUMMARY_RE.search(output)
        if summary is not None:
            if int(summary.group("mode")) != 1:
                raise RuntimeError("unexpected mode, expected cpp mode=1")
            inform = int(summary.group("inform"))
            if inform < 0 or inform > 8:
                raise RuntimeError(f"unexpected inform={inform}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
