#!/usr/bin/env python3
"""Validate the native-AD CSV artifact without redefining scientific tolerances."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--expected-cases", type=int, default=12)
    args = parser.parse_args()

    with args.csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != args.expected_cases:
        raise SystemExit(f"expected {args.expected_cases} AD cases, found {len(rows)}")
    identifiers = {row["case_id"] for row in rows}
    if len(identifiers) != len(rows):
        raise SystemExit("AD benchmark case identifiers are not unique")

    speedups: list[float] = []
    for row in rows:
        absolute = float(row["max_abs_error"])
        relative = float(row["max_rel_error"])
        gap = float(row["spectral_gap"])
        residual = float(row["residual"])
        analytic = int(row["analytic_eigensolves"])
        finite_difference = int(row["finite_difference_eigensolves"])
        parameters = int(row["parameters"])
        speedup = float(row["speedup"])
        if absolute > 2.0e-5 and relative > 2.0e-4:
            raise SystemExit(f"{row['case_id']}: gradient mismatch")
        if not math.isfinite(gap) or gap <= 0.0:
            raise SystemExit(f"{row['case_id']}: unresolved ground-state gap")
        if not math.isfinite(residual) or residual > 1.0e-8:
            raise SystemExit(f"{row['case_id']}: residual {residual} exceeds gate")
        if analytic != 1 or finite_difference != 2 * parameters:
            raise SystemExit(f"{row['case_id']}: incorrect eigensolve accounting")
        if not math.isfinite(speedup) or speedup <= 0.0:
            raise SystemExit(f"{row['case_id']}: invalid timing ratio")
        speedups.append(speedup)

    geometric_mean = math.exp(sum(math.log(value) for value in speedups) / len(speedups))
    print(f"validated {len(rows)} AD workflows; geometric-mean speedup={geometric_mean:.3f}x")
    if geometric_mean <= 1.0:
        raise SystemExit("native AD did not beat central finite differences in aggregate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
