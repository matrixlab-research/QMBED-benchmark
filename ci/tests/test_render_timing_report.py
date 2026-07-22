from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from render_timing_report import render_benchmarks


class TimingReportTests(unittest.TestCase):
    def write_csv(self, directory: Path, language: str, median: float) -> Path:
        path = directory / f"{language}.csv"
        fields = [
            "language",
            "suite",
            "case_id",
            "benchmark",
            "comparison",
            "storage",
            "validation",
            "median_seconds",
            "p25_seconds",
            "p75_seconds",
        ]
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerow(
                {
                    "language": language,
                    "suite": "paper",
                    "case_id": "paper_pxp_revival_l24",
                    "benchmark": "PXP constrained-state revival",
                    "comparison": "current_backend",
                    "storage": "csc",
                    "validation": "physical_invariant",
                    "median_seconds": median,
                    "p25_seconds": median * 0.9,
                    "p75_seconds": median * 1.1,
                }
            )
        return path

    def test_optional_rust_candidate_renders_without_a_julia_csv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            report = render_benchmarks(
                [
                    self.write_csv(directory, "python", 0.2),
                    self.write_csv(directory, "rust", 0.1),
                ]
            )

        self.assertIn("Python QuSpin baseline vs Rust candidate", report)
        self.assertIn("PXP constrained-state revival", report)
        self.assertIn("2.00×", report)
        self.assertNotIn("Missing paper benchmark counterparts", report)


if __name__ == "__main__":
    unittest.main()
