#!/usr/bin/env python3
"""Render correctness, timing, and peak-memory evidence from one hosted runner."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from render_paper_workflow_chart import CASE_ORDER


LANGUAGES = ("python", "julia", "rust")


@dataclass(frozen=True)
class Timing:
    label: str
    seconds: float


def load_timings(path: Path, language: str) -> dict[str, Timing]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    timings = {}
    for row in rows:
        if row.get("suite") != "paper":
            continue
        if row.get("language") != language:
            raise ValueError(f"{path} does not contain {language} timing rows")
        case_id = row["case_id"]
        if case_id in timings:
            raise ValueError(f"{path} contains duplicate case {case_id}")
        if row.get("supported", "true").lower() != "true":
            raise ValueError(f"{language} marks {case_id} unsupported")
        if row.get("validation") != "passed":
            raise ValueError(f"{language} did not validate {case_id}")
        timings[case_id] = Timing(
            row.get("benchmark", case_id),
            float(row["median_seconds"]),
        )
    expected = set(CASE_ORDER)
    if set(timings) != expected:
        missing = sorted(expected - set(timings))
        extra = sorted(set(timings) - expected)
        raise ValueError(
            f"{language} workflow set differs: missing={missing}, extra={extra}"
        )
    return timings


def load_memory(path: Path, language: str) -> tuple[int, float]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f"{path} must contain one memory row")
    row = rows[0]
    if row.get("language") != language:
        raise ValueError(f"{path} does not describe {language}")
    if row.get("metric") != "peak_resident_bytes" or row.get("unit") != "bytes":
        raise ValueError(f"{path} has an unexpected memory metric")
    if int(row.get("exit_code", "1")) != 0:
        raise ValueError(f"{language} benchmark command failed")
    return int(row["value"]), float(row["wall_seconds"])


def render(
    timing_paths: dict[str, Path],
    memory_paths: dict[str, Path],
) -> str:
    timings = {
        language: load_timings(timing_paths[language], language)
        for language in LANGUAGES
    }
    memories = {
        language: load_memory(memory_paths[language], language)
        for language in LANGUAGES
    }
    lines = [
        "# Same-runner Python / Julia / Rust paper benchmark",
        "",
        "Every row passed its language-specific physical invariant before timing. "
        "All processes ran on one hosted runner under the same single-thread "
        "environment; lower is better.",
        "",
        "| Workflow | Python (s) | Julia (s) | Rust (s) | Fastest | Rust / Python |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for case_id in CASE_ORDER:
        values = {
            language: timings[language][case_id].seconds
            for language in LANGUAGES
        }
        label = timings["python"][case_id].label
        fastest = min(values, key=values.get)
        rust_speedup = values["python"] / values["rust"]
        lines.append(
            f"| {label} | {values['python']:.6f} | {values['julia']:.6f} | "
            f"{values['rust']:.6f} | {fastest.title()} | {rust_speedup:.2f}x |"
        )
    lines.extend(
        [
            "",
            "Peak RSS is measured for each complete twelve-workflow process, "
            "so it includes the language runtime and loaded numerical libraries.",
            "",
            "| Language | Peak RSS (MiB) | Suite wall time (s) |",
            "|---|---:|---:|",
        ]
    )
    for language in LANGUAGES:
        peak_bytes, wall_seconds = memories[language]
        lines.append(
            f"| {language.title()} | {peak_bytes / 2**20:.1f} | {wall_seconds:.3f} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python-timing", type=Path, required=True)
    parser.add_argument("--julia-timing", type=Path, required=True)
    parser.add_argument("--rust-timing", type=Path, required=True)
    parser.add_argument("--python-memory", type=Path, required=True)
    parser.add_argument("--julia-memory", type=Path, required=True)
    parser.add_argument("--rust-memory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = render(
        {
            "python": args.python_timing,
            "julia": args.julia_timing,
            "rust": args.rust_timing,
        },
        {
            "python": args.python_memory,
            "julia": args.julia_memory,
            "rust": args.rust_memory,
        },
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
