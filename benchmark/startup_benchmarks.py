#!/usr/bin/env python3
"""Measure package load in fresh processes after both packages are installed."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path


SAMPLES = 9


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def sample(command: list[str]) -> list[float]:
    subprocess.run(
        command,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    values = []
    for _ in range(SAMPLES):
        started = time.perf_counter_ns()
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        values.append((time.perf_counter_ns() - started) / 1.0e9)
    return values


def row(language: str, values: list[float]) -> dict[str, object]:
    mean = statistics.fmean(values)
    return {
        "language": language,
        "benchmark": "package_load_fresh_process",
        "category": "startup",
        "parameters": "precompiled=true;fresh_process=true",
        "samples": len(values),
        "iterations_per_sample": 1,
        "median_seconds": statistics.median(values),
        "mean_seconds": mean,
        "stdev_seconds": statistics.stdev(values),
        "min_seconds": min(values),
        "p05_seconds": percentile(values, 0.05),
        "p25_seconds": percentile(values, 0.25),
        "p75_seconds": percentile(values, 0.75),
        "p95_seconds": percentile(values, 0.95),
        "max_seconds": max(values),
        "runtime": language,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--project", type=Path, default=Path.cwd())
    args = parser.parse_args()
    julia = shutil.which("julia")
    if julia is None:
        raise RuntimeError("julia executable not found")
    rows = [
        row("python", sample([sys.executable, "-c", "import quspin"])),
        row(
            "julia",
            sample(
                [
                    julia,
                    "--startup-file=no",
                    f"--project={args.project.resolve()}",
                    "-e",
                    "using QuSpin",
                ]
            ),
        ),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
