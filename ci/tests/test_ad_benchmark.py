from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HEADER = [
    "case_id",
    "family",
    "dimension",
    "parameters",
    "max_abs_error",
    "max_rel_error",
    "spectral_gap",
    "residual",
    "analytic_seconds",
    "finite_difference_seconds",
    "speedup",
    "analytic_eigensolves",
    "finite_difference_eigensolves",
]


def write_fixture(path: Path, speedup: float = 2.0) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        for index in range(12):
            writer.writerow(
                {
                    "case_id": f"ad_case_{index}",
                    "family": "fixture",
                    "dimension": 32,
                    "parameters": 3,
                    "max_abs_error": "1e-9",
                    "max_rel_error": "1e-9",
                    "spectral_gap": "0.2",
                    "residual": "1e-12",
                    "analytic_seconds": "0.1",
                    "finite_difference_seconds": str(0.1 * speedup),
                    "speedup": str(speedup),
                    "analytic_eigensolves": 1,
                    "finite_difference_eigensolves": 6,
                }
            )


class NativeAdBenchmarkToolingTests(unittest.TestCase):
    def test_validator_accepts_twelve_correct_and_faster_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ad.csv"
            write_fixture(path)
            result = subprocess.run(
                [sys.executable, str(ROOT / "ci/check_ad_benchmark.py"), str(path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("geometric-mean speedup=2.000x", result.stdout)

    def test_validator_rejects_no_aggregate_advantage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ad.csv"
            write_fixture(path, speedup=0.9)
            result = subprocess.run(
                [sys.executable, str(ROOT / "ci/check_ad_benchmark.py"), str(path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("did not beat", result.stderr)

    def test_renderer_emits_all_cases_and_speedup_labels(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            csv_path = Path(directory) / "ad.csv"
            svg_path = Path(directory) / "ad.svg"
            write_fixture(csv_path)
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "ci/render_ad_benchmark.py"),
                    str(csv_path),
                    str(svg_path),
                ],
                check=True,
            )
            svg = svg_path.read_text(encoding="utf-8")
            for index in range(12):
                self.assertIn(f"case {index}", svg)
            self.assertEqual(svg.count("2.00×"), 12)


if __name__ == "__main__":
    unittest.main()
