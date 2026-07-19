#!/usr/bin/env python3
"""Compare oracle JSON structurally while tolerating only BLAS-scale float tails."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


RELATIVE_TOLERANCE = 5.0e-12
ABSOLUTE_TOLERANCE = 5.0e-12


def compare(expected: Any, actual: Any, path: str, errors: list[str]) -> None:
    if isinstance(expected, bool) or isinstance(actual, bool):
        if expected is not actual:
            errors.append(f"{path}: expected {expected!r}, got {actual!r}")
        return

    if isinstance(expected, dict) and isinstance(actual, dict):
        expected_keys = set(expected)
        actual_keys = set(actual)
        for key in sorted(expected_keys - actual_keys):
            errors.append(f"{path}: missing key {key!r}")
        for key in sorted(actual_keys - expected_keys):
            errors.append(f"{path}: unexpected key {key!r}")
        for key in sorted(expected_keys & actual_keys):
            compare(expected[key], actual[key], f"{path}.{key}", errors)
        return

    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            errors.append(
                f"{path}: expected list length {len(expected)}, got {len(actual)}"
            )
            return
        for index, (expected_value, actual_value) in enumerate(
            zip(expected, actual)
        ):
            compare(
                expected_value,
                actual_value,
                f"{path}[{index}]",
                errors,
            )
        return

    if isinstance(expected, int) and isinstance(actual, int):
        if expected != actual:
            errors.append(f"{path}: expected {expected}, got {actual}")
        return

    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        if not math.isclose(
            float(expected),
            float(actual),
            rel_tol=RELATIVE_TOLERANCE,
            abs_tol=ABSOLUTE_TOLERANCE,
        ):
            errors.append(f"{path}: expected {expected:.17g}, got {actual:.17g}")
        return

    if type(expected) is not type(actual) or expected != actual:
        errors.append(f"{path}: expected {expected!r}, got {actual!r}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("expected", type=Path)
    parser.add_argument("actual", type=Path)
    args = parser.parse_args()
    expected = json.loads(args.expected.read_text())
    actual = json.loads(args.actual.read_text())
    errors: list[str] = []
    compare(expected, actual, "$", errors)
    if errors:
        details = "\n".join(errors[:50])
        remaining = len(errors) - 50
        if remaining > 0:
            details += f"\n... and {remaining} more differences"
        raise SystemExit(
            "oracle comparison failed "
            f"(rtol={RELATIVE_TOLERANCE:g}, atol={ABSOLUTE_TOLERANCE:g}):\n"
            + details
        )
    print(
        "oracle comparison passed "
        f"(rtol={RELATIVE_TOLERANCE:g}, atol={ABSOLUTE_TOLERANCE:g})"
    )


if __name__ == "__main__":
    main()
