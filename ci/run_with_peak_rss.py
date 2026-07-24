#!/usr/bin/env python3
"""Run one benchmark command and record its wall time and peak resident set."""

from __future__ import annotations

import argparse
import csv
import platform
import resource
import shlex
import subprocess
import sys
import time
from pathlib import Path


FIELDS = (
    "language",
    "suite",
    "metric",
    "value",
    "unit",
    "wall_seconds",
    "exit_code",
    "command",
)


def peak_resident_bytes() -> int:
    maximum = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return int(maximum if platform.system() == "Darwin" else maximum * 1024)


def run(
    command: list[str],
    *,
    language: str,
    output: Path,
    stdout: Path | None = None,
) -> int:
    if not command:
        raise ValueError("a benchmark command is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    if stdout is not None:
        stdout.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    if stdout is None:
        completed = subprocess.run(command, check=False)
    else:
        with stdout.open("w") as stream:
            completed = subprocess.run(
                command,
                check=False,
                stdout=stream,
                text=True,
            )
    wall_seconds = time.perf_counter() - started
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerow(
            {
                "language": language,
                "suite": "paper",
                "metric": "peak_resident_bytes",
                "value": peak_resident_bytes(),
                "unit": "bytes",
                "wall_seconds": f"{wall_seconds:.12f}",
                "exit_code": completed.returncode,
                "command": shlex.join(command),
            }
        )
    return completed.returncode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--language", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stdout", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    raise SystemExit(
        run(
            command,
            language=args.language,
            output=args.output,
            stdout=args.stdout,
        )
    )


if __name__ == "__main__":
    main()
