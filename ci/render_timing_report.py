#!/usr/bin/env python3
"""Render test-suite and cross-language timing CSV files as Markdown."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def milliseconds(value: str) -> str:
    return f"{float(value) * 1_000:.3f}"


def mib(value: str) -> str:
    return f"{int(value) / (1024 * 1024):.2f}"


def render_tests(path: Path) -> str:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    suite = next(row for row in rows if row["path"] == "__suite__")
    files = [row for row in rows if row["path"] != "__suite__"]
    files.sort(key=lambda row: float(row["wall_seconds"]), reverse=True)
    lines = [
        "## Private verification timing",
        "",
        (
            f"Total wall time: **{float(suite['wall_seconds']):.3f} s**; "
            f"reported compile time: **{float(suite['compile_seconds']):.3f} s**; "
            f"allocations: **{mib(suite['allocated_bytes'])} MiB**."
        ),
        "",
        "| Test file | Scope | Wall (ms) | Compile (ms) | GC (ms) | Allocated (MiB) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in files:
        lines.append(
            "| `{path}` | {scope} | {wall} | {compile} | {gc} | {allocated} |".format(
                path=row["path"],
                scope=row["scope"],
                wall=milliseconds(row["wall_seconds"]),
                compile=milliseconds(row["compile_seconds"]),
                gc=milliseconds(row["gc_seconds"]),
                allocated=mib(row["allocated_bytes"]),
            )
        )
    return "\n".join(lines) + "\n"


def read_benchmarks(paths: list[Path]) -> dict[str, dict[str, dict[str, str]]]:
    cases: dict[str, dict[str, dict[str, str]]] = {}
    for path in paths:
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                cases.setdefault(row["benchmark"], {})[row["language"]] = row
    return cases


def render_benchmarks(paths: list[Path]) -> str:
    cases = read_benchmarks(paths)
    lines = [
        "## Python QuSpin baseline vs Julia candidate",
        "",
        (
            "Warm-operation timings exclude JIT warm-up. Both implementations run "
            "on the same GitHub Actions runner with one Julia and BLAS thread. "
            "Speedup is `Python median / Julia median`; values above 1 mean Julia is faster."
        ),
        "",
        (
            "`controlled` rows compare the same dense or CSC representation. "
            "`current_backend` rows intentionally retain each package's present "
            "public-API storage choice."
        ),
        "",
        "| Workload | Mode | Python storage | Julia storage | Python median (ms) | Julia median (ms) | Julia IQR (ms) | Speedup | Julia allocation |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    notes = []
    for name, languages in cases.items():
        if set(languages) != {"python", "julia"}:
            continue
        python = languages["python"]
        julia = languages["julia"]
        julia_supported = julia.get("supported", "true").lower() == "true"
        if julia_supported:
            speedup = (
                float(python["median_seconds"]) / float(julia["median_seconds"])
            )
            julia_ms = milliseconds(julia["median_seconds"])
            julia_iqr = (
                f"{milliseconds(julia['p25_seconds'])}–"
                f"{milliseconds(julia['p75_seconds'])}"
            )
            speedup_text = f"{speedup:.2f}×"
            allocation = julia.get("median_allocated_bytes")
            allocation_text = (
                f"{int(allocation) / 1024:.1f} KiB" if allocation else "n/a"
            )
        else:
            julia_ms = "unsupported"
            julia_iqr = "n/a"
            speedup_text = "n/a"
            allocation_text = "n/a"
            if julia.get("note"):
                notes.append(f"- `{name}`: {julia['note']}")
        lines.append(
            "| `{name}` | {comparison} | {python_storage} | {julia_storage} | "
            "{python_ms} | {julia_ms} | {julia_iqr} | {speedup} | {allocation} |".format(
                name=name,
                comparison=julia.get("comparison", "unspecified"),
                python_storage=python.get("storage", "unspecified"),
                julia_storage=julia.get("storage", "unspecified"),
                python_ms=milliseconds(python["median_seconds"]),
                julia_ms=julia_ms,
                julia_iqr=julia_iqr,
                speedup=speedup_text,
                allocation=allocation_text,
            )
        )
    if notes:
        lines.extend(["", "Capability gaps:", "", *notes])
    lines.extend(
        [
            "",
            "Operation benchmarks use 15 samples; fresh-process loading uses 9. "
            "The uploaded CSV files retain aggregate statistics for both: "
            "minimum, p05, p25, median, mean, p75, p95, maximum, standard deviation, "
            "and iterations per sample.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    tests = subparsers.add_parser("tests")
    tests.add_argument("path", type=Path)
    tests.add_argument("--output", type=Path)
    benchmarks = subparsers.add_parser("benchmarks")
    benchmarks.add_argument("paths", type=Path, nargs="+")
    benchmarks.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = (
        render_tests(args.path)
        if args.command == "tests"
        else render_benchmarks(args.paths)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
