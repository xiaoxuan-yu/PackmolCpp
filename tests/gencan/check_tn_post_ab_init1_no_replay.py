#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import tempfile


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    probe = repo / "build" / "tests" / "test_gencan_ab_probe"
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
                "PACKMOL_GENCAN_IMPL": "ab",
                "PACKMOL_GENCAN_NUMERIC_CPP": "0",
                "PACKMOL_GENCAN_AB_COMPARE": "1",
                "PACKMOL_GENCAN_DEBUG": "1",
                "PACKMOL_GENCAN_TN_POST_RETRY_SPG": "0",
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
        if completed.returncode != 0:
            raise RuntimeError(f"probe failed\n{output}")

        if "[gencan-cpp-tn-post-cpp-replay] depth=" in output:
            raise RuntimeError(
                "tn_post AB init1 path unexpectedly entered recursive cpp replay\n"
                f"{output}"
            )
        if "[gencan-cpp-tn-post-cpp-replay-gate]" not in output:
            raise RuntimeError(
                "tn_post AB init1 path did not emit replay gate diagnostics\n"
                f"{output}"
            )
        if "replay_enabled=0" not in output:
            raise RuntimeError(
                "tn_post AB init1 replay gate was expected to be disabled\n"
                f"{output}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
