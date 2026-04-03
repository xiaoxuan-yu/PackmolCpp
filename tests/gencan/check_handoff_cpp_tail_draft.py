#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
import tempfile


PROBES = (
    "test_gencan_ab_probe",
    "test_gencan_spg_post_ab_probe",
)

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
        env.update({
            "PACKMOL_GENCAN_IMPL": "cpp",
            "PACKMOL_GENCAN_NUMERIC_CPP": "0",
            "PACKMOL_GENCAN_DEBUG": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_UNSAFE": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_MAX_DEPTH": "1",
            "PACKMOL_GENCAN_TN_POST_CPP_HANDOFF_CPP_REPLAY": "1",
            "PACKMOL_GENCAN_SPG_POST_CPP_HANDOFF": "1",
            "PACKMOL_GENCAN_HANDOFF_CPP_TAIL": "1",
            "PACKMOL_GENCAN_HANDOFF_CPP_TAIL_ALLOW_TN_POST": "1",
            "PACKMOL_GENCAN_HANDOFF_CPP_TAIL_ALLOW_SPG_POST": "1",
        })

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
            if "[gencan-cpp-fortran-tail] reason=spg_post_nonterminal" in output:
                raise RuntimeError(f"{probe_name}: spg_post_nonterminal unexpectedly used fortran-tail")
            if probe_name == "test_gencan_spg_post_ab_probe" and "[gencan-cpp-handoff-cpp-tail]" not in output:
                raise RuntimeError(f"{probe_name}: missing handoff-cpp-tail marker")
            if SUMMARY_RE.search(output) is None:
                raise RuntimeError(f"{probe_name}: missing summary line")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
