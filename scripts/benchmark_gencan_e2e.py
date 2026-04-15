#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import re
import statistics
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass
from typing import Any


DEFAULT_CASES = [
    "water_box.inp",
    "water_box_pbc.inp",
    "bilayer.inp",
    "mixture.inp",
]
MODES = ("fortran", "cpp")


@dataclass
class ArtifactRecord:
    relative_path: str
    size: int
    sha256: str


@dataclass
class RunSummary:
    returncode: int
    wall_time_s: float
    loop_count: int
    success: bool
    final_objective: float | None
    final_target_violation: float | None
    final_constraint_violation: float | None
    pair_penalty_sum: float | None
    constraint_penalty_sum: float | None
    pair_penalty_count: int | None
    constraint_penalty_count: int | None
    normalized_log_sha256: str
    artifacts: list[ArtifactRecord]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark and compare Packmol GENCAN fortran/cpp end-to-end behavior."
    )
    parser.add_argument(
        "--build-dir",
        default="build",
        help="Build directory that contains the packmol executable. Default: %(default)s",
    )
    parser.add_argument(
        "--binary",
        default="packmol",
        help="Executable name inside the build directory. Default: %(default)s",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=DEFAULT_CASES,
        help="Input files under testing/input_files/. Default: %(default)s",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=1,
        help="Warmup runs per case/mode before measurement. Default: %(default)s",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Measured repeats per case/mode. Default: %(default)s",
    )
    parser.add_argument(
        "--json-out",
        help="Optional path to write the full JSON report.",
    )
    return parser.parse_args()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_log(text: str) -> str:
    lines: list[str] = []
    for raw_line in text.splitlines():
        line = re.sub(
            r"Running time:\s+[-+0-9.Ee]+\s+seconds\.",
            "Running time: <normalized> seconds.",
            raw_line,
        )
        line = re.sub(
            r"Starting GENCAN loop:\s+\d+",
            "Starting GENCAN loop: <n>",
            line,
        )
        lines.append(line.rstrip())
    return "\n".join(lines) + "\n"


def collect_artifacts(workdir: pathlib.Path) -> list[ArtifactRecord]:
    artifacts: list[ArtifactRecord] = []
    for root, dirs, files in os.walk(workdir, topdown=True, followlinks=False):
        root_path = pathlib.Path(root)
        dirs[:] = [name for name in dirs if name not in {"input_files", "structure_files"}]
        for name in sorted(files):
            if name == "run.log":
                continue
            path = root_path / name
            relative_path = path.relative_to(workdir).as_posix()
            artifacts.append(
                ArtifactRecord(
                    relative_path=relative_path,
                    size=path.stat().st_size,
                    sha256=hash_file(path),
                )
            )
    artifacts.sort(key=lambda item: item.relative_path)
    return artifacts


def parse_last_float(pattern: str, text: str) -> float | None:
    matches = re.findall(pattern, text)
    if not matches:
        return None
    return float(matches[-1])


def parse_last_int(pattern: str, text: str) -> int | None:
    matches = re.findall(pattern, text)
    if not matches:
        return None
    return int(matches[-1])


def parse_run_summary(workdir: pathlib.Path, returncode: int, wall_time_s: float) -> RunSummary:
    log_text = (workdir / "run.log").read_text(encoding="utf-8", errors="replace")
    normalized_log = normalize_log(log_text)
    return RunSummary(
        returncode=returncode,
        wall_time_s=wall_time_s,
        loop_count=log_text.count("Starting GENCAN loop"),
        success="Success!" in log_text,
        final_objective=parse_last_float(r"Final objective function value:\s*([-+.\dEe]+)", log_text),
        final_target_violation=parse_last_float(
            r"Maximum violation of target distance:\s*([-+.\dEe]+)",
            log_text,
        ),
        final_constraint_violation=parse_last_float(
            r"Maximum violation of the constraints:\s*([-+.\dEe]+)",
            log_text,
        ),
        pair_penalty_sum=parse_last_float(r"Pair penalty sum:\s*([-+.\dEe]+)", log_text),
        constraint_penalty_sum=parse_last_float(r"Constraint penalty sum:\s*([-+.\dEe]+)", log_text),
        pair_penalty_count=parse_last_int(r"Active pair penalty count:\s*(\d+)", log_text),
        constraint_penalty_count=parse_last_int(r"Active constraint penalty count:\s*(\d+)", log_text),
        normalized_log_sha256=sha256_bytes(normalized_log.encode("utf-8")),
        artifacts=collect_artifacts(workdir),
    )


def run_once(
    binary_path: pathlib.Path,
    input_root: pathlib.Path,
    structure_root: pathlib.Path,
    case: str,
    mode: str,
) -> RunSummary:
    with tempfile.TemporaryDirectory() as tempdir:
        workdir = pathlib.Path(tempdir)
        (workdir / "input_files").symlink_to(input_root)
        (workdir / "structure_files").symlink_to(structure_root)
        env = os.environ.copy()
        env["PACKMOL_GENCAN_IMPL"] = mode
        env["PACKMOL_REPORT_PENALTY_STATS"] = "1"
        input_path = workdir / "input_files" / case
        start = time.perf_counter()
        with input_path.open("rb") as stdin_handle, (workdir / "run.log").open("wb") as log_handle:
            process = subprocess.run(
                [str(binary_path)],
                stdin=stdin_handle,
                cwd=workdir,
                env=env,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=False,
            )
        wall_time_s = time.perf_counter() - start
        return parse_run_summary(workdir, process.returncode, wall_time_s)


def compare_behavior(fortran_run: RunSummary, cpp_run: RunSummary) -> dict[str, Any]:
    fortran_artifacts = {item.relative_path: item for item in fortran_run.artifacts}
    cpp_artifacts = {item.relative_path: item for item in cpp_run.artifacts}
    artifact_paths_match = sorted(fortran_artifacts) == sorted(cpp_artifacts)
    differing_artifacts: list[dict[str, Any]] = []
    for relative_path in sorted(set(fortran_artifacts) & set(cpp_artifacts)):
        left = fortran_artifacts[relative_path]
        right = cpp_artifacts[relative_path]
        if left.sha256 != right.sha256 or left.size != right.size:
            differing_artifacts.append(
                {
                    "relative_path": relative_path,
                    "fortran_size": left.size,
                    "cpp_size": right.size,
                    "fortran_sha256": left.sha256,
                    "cpp_sha256": right.sha256,
                }
            )

    return {
        "same_returncode": fortran_run.returncode == cpp_run.returncode,
        "same_success_flag": fortran_run.success == cpp_run.success,
        "same_loop_count": fortran_run.loop_count == cpp_run.loop_count,
        "same_artifact_paths": artifact_paths_match,
        "same_artifact_bytes": artifact_paths_match and not differing_artifacts,
        "same_normalized_log": fortran_run.normalized_log_sha256 == cpp_run.normalized_log_sha256,
        "differing_artifacts": differing_artifacts,
        "fortran_final_objective": fortran_run.final_objective,
        "cpp_final_objective": cpp_run.final_objective,
        "fortran_final_target_violation": fortran_run.final_target_violation,
        "cpp_final_target_violation": cpp_run.final_target_violation,
        "fortran_final_constraint_violation": fortran_run.final_constraint_violation,
        "cpp_final_constraint_violation": cpp_run.final_constraint_violation,
        "fortran_pair_penalty_sum": fortran_run.pair_penalty_sum,
        "cpp_pair_penalty_sum": cpp_run.pair_penalty_sum,
        "fortran_constraint_penalty_sum": fortran_run.constraint_penalty_sum,
        "cpp_constraint_penalty_sum": cpp_run.constraint_penalty_sum,
        "fortran_pair_penalty_count": fortran_run.pair_penalty_count,
        "cpp_pair_penalty_count": cpp_run.pair_penalty_count,
        "fortran_constraint_penalty_count": fortran_run.constraint_penalty_count,
        "cpp_constraint_penalty_count": cpp_run.constraint_penalty_count,
    }


def benchmark_case(
    binary_path: pathlib.Path,
    input_root: pathlib.Path,
    structure_root: pathlib.Path,
    case: str,
    warmup: int,
    repeats: int,
) -> dict[str, Any]:
    case_report: dict[str, Any] = {"modes": {}, "behavior": None}
    baseline_runs: dict[str, RunSummary] = {}

    for mode in MODES:
        for _ in range(warmup):
            warmup_run = run_once(binary_path, input_root, structure_root, case, mode)
            if warmup_run.returncode != 0:
                raise RuntimeError(f"{case} warmup failed for {mode} with return code {warmup_run.returncode}")

        measured_runs: list[RunSummary] = []
        for _ in range(repeats):
            run = run_once(binary_path, input_root, structure_root, case, mode)
            if run.returncode != 0:
                raise RuntimeError(f"{case} benchmark failed for {mode} with return code {run.returncode}")
            measured_runs.append(run)

        baseline_runs[mode] = measured_runs[0]
        wall_times = [run.wall_time_s for run in measured_runs]
        loop_counts = [run.loop_count for run in measured_runs]
        case_report["modes"][mode] = {
            "median_wall_time_s": statistics.median(wall_times),
            "min_wall_time_s": min(wall_times),
            "max_wall_time_s": max(wall_times),
            "wall_time_runs_s": wall_times,
            "loop_counts": loop_counts,
            "baseline_run": asdict(measured_runs[0]),
        }

    case_report["behavior"] = compare_behavior(baseline_runs["fortran"], baseline_runs["cpp"])
    return case_report


def print_human_report(report: dict[str, Any]) -> None:
    print("GENCAN E2E Benchmark")
    print()
    print(
        f"binary={report['binary']} warmup={report['warmup']} repeats={report['repeats']} "
        f"cases={','.join(report['cases'].keys())}"
    )
    print()
    print(
        f"{'case':20} {'fortran(s)':>12} {'cpp(s)':>12} {'cpp/fortran':>12} "
        f"{'loops f/cpp':>12} {'behavior':>12}"
    )
    for case, case_report in report["cases"].items():
        fortran_info = case_report["modes"]["fortran"]
        cpp_info = case_report["modes"]["cpp"]
        behavior = case_report["behavior"]
        ratio = cpp_info["median_wall_time_s"] / fortran_info["median_wall_time_s"]
        loops_text = f"{fortran_info['baseline_run']['loop_count']}/{cpp_info['baseline_run']['loop_count']}"
        behavior_label = (
            "byte-equal"
            if behavior["same_artifact_bytes"] and behavior["same_normalized_log"]
            else "same-status"
            if behavior["same_returncode"] and behavior["same_success_flag"] and behavior["same_artifact_paths"]
            else "changed"
        )
        print(
            f"{case:20} {fortran_info['median_wall_time_s']:12.3f} {cpp_info['median_wall_time_s']:12.3f} "
            f"{ratio:12.3f} {loops_text:>12} "
            f"{behavior_label:>12}"
        )
    print()
    overall = report["overall"]
    print(
        f"overall_fortran={overall['fortran_median_total_s']:.3f}s "
        f"overall_cpp={overall['cpp_median_total_s']:.3f}s "
        f"ratio={overall['cpp_over_fortran']:.3f}"
    )
    print()
    print("Behavior details")
    for case, case_report in report["cases"].items():
        behavior = case_report["behavior"]
        print(
            f"- {case}: returncode={behavior['same_returncode']} success={behavior['same_success_flag']} "
            f"loops={behavior['same_loop_count']} artifact_paths={behavior['same_artifact_paths']} "
            f"artifact_bytes={behavior['same_artifact_bytes']} normalized_log={behavior['same_normalized_log']} "
            f"final_obj={behavior['fortran_final_objective']}/{behavior['cpp_final_objective']} "
            f"final_constraint={behavior['fortran_final_constraint_violation']}/{behavior['cpp_final_constraint_violation']} "
            f"pair_sum={behavior['fortran_pair_penalty_sum']}/{behavior['cpp_pair_penalty_sum']} "
            f"constraint_sum={behavior['fortran_constraint_penalty_sum']}/{behavior['cpp_constraint_penalty_sum']} "
            f"pair_count={behavior['fortran_pair_penalty_count']}/{behavior['cpp_pair_penalty_count']} "
            f"constraint_count={behavior['fortran_constraint_penalty_count']}/{behavior['cpp_constraint_penalty_count']}"
        )


def main() -> int:
    args = parse_args()
    repo_root = pathlib.Path(__file__).resolve().parent.parent
    binary_path = repo_root / args.build_dir / args.binary
    input_root = repo_root / "testing" / "input_files"
    structure_root = repo_root / "testing" / "structure_files"

    if not binary_path.exists():
        raise SystemExit(f"packmol binary not found: {binary_path}")

    report: dict[str, Any] = {
        "binary": str(binary_path),
        "warmup": args.warmup,
        "repeats": args.repeats,
        "cases": {},
        "overall": {},
    }

    for case in args.cases:
        report["cases"][case] = benchmark_case(
            binary_path=binary_path,
            input_root=input_root,
            structure_root=structure_root,
            case=case,
            warmup=args.warmup,
            repeats=args.repeats,
        )

    fortran_total = sum(
        case_report["modes"]["fortran"]["median_wall_time_s"]
        for case_report in report["cases"].values()
    )
    cpp_total = sum(
        case_report["modes"]["cpp"]["median_wall_time_s"]
        for case_report in report["cases"].values()
    )
    report["overall"] = {
        "fortran_median_total_s": fortran_total,
        "cpp_median_total_s": cpp_total,
        "cpp_over_fortran": cpp_total / fortran_total,
    }

    if args.json_out:
        json_path = pathlib.Path(args.json_out)
        json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print_human_report(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
