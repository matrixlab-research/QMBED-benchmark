from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from render_paper_workflow_chart import CASE_ORDER, load_pairs, render_svg


class PaperWorkflowChartTests(unittest.TestCase):
    def write_csv(
        self,
        directory: Path,
        language: str,
        *,
        case_ids: tuple[str, ...] = CASE_ORDER,
    ) -> Path:
        path = directory / f"{language}.csv"
        fields = [
            "language",
            "suite",
            "case_id",
            "benchmark",
            "supported",
            "median_seconds",
        ]
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for index, case_id in enumerate(case_ids, start=1):
                writer.writerow(
                    {
                        "language": language,
                        "suite": "paper",
                        "case_id": case_id,
                        "benchmark": f"Workflow {index}",
                        "supported": "true",
                        "median_seconds": 0.01 * index * (2 if language == "python" else 1),
                    }
                )
        return path

    def test_renders_twelve_paired_bars_with_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            pairs = load_pairs(
                self.write_csv(directory, "python"),
                self.write_csv(directory, "julia"),
            )
            svg = render_svg(
                pairs,
                run_id="42",
                candidate_repository="matrixlab-research/QuSpin.jl",
                candidate_ref="0123456789abcdef",
            )

        self.assertEqual(svg.count('class="python-bar"'), 12)
        self.assertEqual(svg.count('class="julia-bar"'), 12)
        self.assertIn("GitHub Actions run 42", svg)
        self.assertIn("0123456789ab", svg)
        self.assertIn("Workflow 12", svg)
        self.assertIn("240.000", svg)
        self.assertIn("120.000", svg)

    def test_rejects_missing_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            python_csv = self.write_csv(directory, "python", case_ids=CASE_ORDER[:-1])
            julia_csv = self.write_csv(directory, "julia")
            with self.assertRaisesRegex(
                ValueError, "missing paper_particle_addition_6x3"
            ):
                load_pairs(python_csv, julia_csv)


if __name__ == "__main__":
    unittest.main()
