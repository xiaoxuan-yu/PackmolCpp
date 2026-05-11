#!/usr/bin/env python3

import re
import subprocess
import sys
from pathlib import Path


def parse_line(pattern: str, text: str) -> re.Match[str]:
    match = re.search(pattern, text)
    if match is None:
        raise SystemExit(f"missing pattern: {pattern}")
    return match


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_cg_direct_compare.py <probe> <input>")

    probe = Path(sys.argv[1]).resolve()
    input_path = Path(sys.argv[2]).resolve()
    workdir = input_path.parent.parent

    proc = subprocess.run(
        [str(probe), str(input_path.relative_to(workdir))],
        cwd=workdir,
        capture_output=True,
        text=True,
        check=False,
    )

    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(proc.returncode)

    output = proc.stdout + proc.stderr

    used_cpp = int(parse_line(r"used_cpp\s+(\d+)", output).group(1))
    if used_cpp != 1:
        raise SystemExit(f"expected used_cpp=1, got {used_cpp}")

    fortran = parse_line(
        r"fortran_inform\s+(\d+)\s+fortran_iter\s+(\d+)\s+fortran_rbdtype\s+(\d+)\s+fortran_rbdind\s+(\d+)",
        output,
    )
    cpp = parse_line(
        r"cpp_inform\s+(\d+)\s+cpp_iter\s+(\d+)\s+cpp_rbdtype\s+(\d+)\s+cpp_rbdind\s+(\d+)",
        output,
    )
    diff = float(parse_line(r"diff_s_max\s+([0-9Ee+\-.]+)", output).group(1))

    if fortran.groups() != cpp.groups():
        raise SystemExit(
            "CG direct compare mismatch: "
            f"fortran={fortran.groups()} cpp={cpp.groups()}"
        )

    if diff > 1.0e-10:
        raise SystemExit(f"diff_s_max too large: {diff}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
