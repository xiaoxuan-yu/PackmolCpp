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
    "test_gencan_entry_stop_ab_probe",
    "test_gencan_entry_stop_gtype2_ab_probe",
    "test_gencan_entry_stop_inform2_ab_probe",
    "test_gencan_entry_stop_inform3_ab_probe",
    "test_gencan_runtime_driver_probe",
)

TAIL_RE = re.compile(r"\[gencan-cpp-fortran-tail\]\s+reason=(?P<reason>[^\s]+)\s+mode=")


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
                "PACKMOL_GENCAN_DEBUG": "1",
            }
        )

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

            for match in TAIL_RE.finditer(output):
                reason = match.group("reason")
                if reason in {"tn_post_nonterminal", "tn_no_free_variables"}:
                    raise RuntimeError(
                        f"{probe_name}: unexpected fortran-tail reason '{reason}' in cpp mode"
                    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
