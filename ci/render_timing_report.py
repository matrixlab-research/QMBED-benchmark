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


def read_benchmarks(
    paths: list[Path],
) -> dict[tuple[str, str], dict[str, dict[str, str]]]:
    cases: dict[tuple[str, str], dict[str, dict[str, str]]] = {}
    for path in paths:
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                suite = row.get("suite") or "micro"
                case_id = row.get("case_id") or row["benchmark"]
                cases.setdefault((suite, case_id), {})[row["language"]] = row
    return cases


def render_benchmarks(paths: list[Path]) -> str:
    cases = read_benchmarks(paths)
    lines = []
    if any("julia" in languages for languages in cases.values()):
        lines.extend(
            [
                "## Python QuSpin baseline vs Julia candidate",
                "",
                (
                    "Warm-operation timings exclude JIT warm-up. Both implementations run "
                    "on the same runner with one Julia and BLAS thread. "
                    "Speedup is `Python median / Julia median`; values above 1 mean Julia is faster."
                ),
                "",
                (
                    "`controlled` rows compare the same dense or named sparse representation. "
                    "`current_backend` rows intentionally retain each package's present "
                    "public-API storage choice."
                ),
                "",
                (
                    "Validation is outside the timed region: `passed` means an explicit "
                    "dimension/shape, matrix-fingerprint, or physical-invariant check passed; "
                    "`smoke` means the warm preflight completed without an exception."
                ),
                "",
                "| Suite | Workload | Mode | Validation | Python storage | Julia storage | Python median (ms) | Julia median (ms) | Julia IQR (ms) | Speedup | Julia allocation |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
    notes = []
    julia_only = []
    rust_pairs = []
    missing_pairs = []
    for (suite, case_id), languages in cases.items():
        has_julia_pair = "python" in languages and "julia" in languages
        has_rust_pair = "python" in languages and "rust" in languages
        if has_rust_pair:
            rust_pairs.append(
                (suite, case_id, languages["python"], languages["rust"])
            )
        if set(languages) == {"julia"}:
            julia_only.append((suite, case_id, languages["julia"]))
        if not has_julia_pair:
            if suite == "paper" and not has_rust_pair:
                missing_pairs.append((suite, case_id, sorted(languages)))
            continue
        python = languages["python"]
        julia = languages["julia"]
        name = julia.get("benchmark") or python.get("benchmark") or case_id
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
        validation = (
            julia.get("validation")
            or python.get("validation")
            or "not_recorded"
        )
        lines.append(
            "| {suite} | `{name}` | {comparison} | {validation} | "
            "{python_storage} | {julia_storage} | "
            "{python_ms} | {julia_ms} | {julia_iqr} | {speedup} | {allocation} |".format(
                suite=suite,
                name=name,
                comparison=julia.get("comparison", "unspecified"),
                validation=validation,
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
    if julia_only:
        lines.extend(
            [
                "",
                "## Julia workflow coverage timings",
                "",
                (
                    "These rows are end-to-end coverage/regression timings, not "
                    "cross-language speedup claims."
                ),
                "",
                "| Suite | Family | Workflow | Parameters | Validation | Median (ms) | IQR (ms) | Allocation |",
                "|---|---:|---|---|---:|---:|---:|---:|",
            ]
        )
        for suite, case_id, julia in julia_only:
            allocation = julia.get("median_allocated_bytes")
            allocation_text = (
                f"{int(allocation) / 1024:.1f} KiB"
                if allocation not in (None, "")
                else "n/a"
            )
            lines.append(
                "| {suite} | {family} | `{name}` | {parameters} | "
                "{validation} | {median} | {p25}–{p75} | {allocation} |".format(
                    suite=suite,
                    family=julia.get("family_id", "n/a"),
                    name=julia.get("benchmark", case_id),
                    parameters=julia.get("parameters", ""),
                    validation=julia.get("validation", "not_recorded"),
                    median=milliseconds(julia["median_seconds"]),
                    p25=milliseconds(julia["p25_seconds"]),
                    p75=milliseconds(julia["p75_seconds"]),
                    allocation=allocation_text,
                )
            )
    if rust_pairs:
        lines.extend(
            [
                "",
                "## Python QuSpin baseline vs Rust candidate",
                "",
                (
                    "The Rust adapter uses the same paper case IDs, physical "
                    "preflights, warm-up count, sample count, and one-thread policy. "
                    "Speedup is `Python median / Rust median`."
                ),
                "",
                "| Suite | Workload | Mode | Validation | Python storage | Rust storage | Python median (ms) | Rust median (ms) | Rust IQR (ms) | Speedup | Rust allocation |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for suite, case_id, python, rust in rust_pairs:
            name = rust.get("benchmark") or python.get("benchmark") or case_id
            rust_supported = rust.get("supported", "true").lower() == "true"
            if rust_supported:
                speedup = (
                    float(python["median_seconds"])
                    / float(rust["median_seconds"])
                )
                rust_ms = milliseconds(rust["median_seconds"])
                rust_iqr = (
                    f"{milliseconds(rust['p25_seconds'])}–"
                    f"{milliseconds(rust['p75_seconds'])}"
                )
                speedup_text = f"{speedup:.2f}×"
                allocation = rust.get("median_allocated_bytes")
                allocation_text = (
                    f"{int(allocation) / 1024:.1f} KiB"
                    if allocation
                    else "n/a"
                )
            else:
                rust_ms = "unsupported"
                rust_iqr = "n/a"
                speedup_text = "n/a"
                allocation_text = "n/a"
            validation = (
                rust.get("validation")
                or python.get("validation")
                or "not_recorded"
            )
            lines.append(
                "| {suite} | `{name}` | {comparison} | {validation} | "
                "{python_storage} | {rust_storage} | "
                "{python_ms} | {rust_ms} | {rust_iqr} | {speedup} | {allocation} |".format(
                    suite=suite,
                    name=name,
                    comparison=rust.get("comparison", "unspecified"),
                    validation=validation,
                    python_storage=python.get("storage", "unspecified"),
                    rust_storage=rust.get("storage", "unspecified"),
                    python_ms=milliseconds(python["median_seconds"]),
                    rust_ms=rust_ms,
                    rust_iqr=rust_iqr,
                    speedup=speedup_text,
                    allocation=allocation_text,
                )
            )
    if missing_pairs:
        lines.extend(["", "## Missing paper benchmark counterparts", ""])
        for suite, case_id, languages in missing_pairs:
            lines.append(
                f"- `{suite}/{case_id}` only has: {', '.join(languages)}."
            )
    lines.extend(
        [
            "",
            "The uploaded CSV files retain aggregate statistics and workflow "
            "benchmarks additionally retain raw samples: "
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
