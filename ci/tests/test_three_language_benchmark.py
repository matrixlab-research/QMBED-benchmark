from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from render_paper_workflow_chart import CASE_ORDER
from render_three_language_report import render
from run_with_peak_rss import run


class ThreeLanguageBenchmarkTests(unittest.TestCase):
    def write_timing(self, directory: Path, language: str, scale: float) -> Path:
        path = directory / f"{language}-timing.csv"
        fields = [
            "language",
            "suite",
            "case_id",
            "benchmark",
            "supported",
            "validation",
            "median_seconds",
        ]
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for index, case_id in enumerate(CASE_ORDER, start=1):
                writer.writerow(
                    {
                        "language": language,
                        "suite": "paper",
                        "case_id": case_id,
                        "benchmark": f"Workflow {index}",
                        "supported": "true",
                        "validation": "passed",
                        "median_seconds": scale * index,
                    }
                )
        return path

    def write_memory(self, directory: Path, language: str) -> Path:
        path = directory / f"{language}-memory.csv"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=[
                    "language",
                    "suite",
                    "metric",
                    "value",
                    "unit",
                    "wall_seconds",
                    "exit_code",
                    "command",
                ],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "language": language,
                    "suite": "paper",
                    "metric": "peak_resident_bytes",
                    "value": 64 * 2**20,
                    "unit": "bytes",
                    "wall_seconds": 2.5,
                    "exit_code": 0,
                    "command": "benchmark",
                }
            )
        return path

    def test_report_requires_validated_twelve_case_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            timing_paths = {
                "python": self.write_timing(directory, "python", 0.03),
                "julia": self.write_timing(directory, "julia", 0.02),
                "rust": self.write_timing(directory, "rust", 0.01),
            }
            memory_paths = {
                language: self.write_memory(directory, language)
                for language in timing_paths
            }
            report = render(timing_paths, memory_paths)
        self.assertIn("Every row passed", report)
        self.assertIn("| Workflow 12 |", report)
        self.assertIn("| Rust | 64.0 | 2.500 |", report)
        self.assertIn("| Rust | 3.00x |", report)

    def test_peak_rss_wrapper_records_a_successful_child(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "memory.csv"
            exit_code = run(
                [sys.executable, "-c", "values = [0] * 10000"],
                language="python",
                output=output,
            )
            with output.open(newline="") as stream:
                row = next(csv.DictReader(stream))
        self.assertEqual(exit_code, 0)
        self.assertEqual(row["metric"], "peak_resident_bytes")
        self.assertGreater(int(row["value"]), 0)
