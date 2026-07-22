from __future__ import annotations

import csv
import re
import sys
import tempfile
import unittest
from pathlib import Path


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from render_paper_workflow_chart import (
    CASE_ORDER,
    _speedup_annotation,
    load_pairs,
    render_svg,
)


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

    def test_renders_twelve_overlaid_pairs_with_speedups_and_provenance(self) -> None:
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
        self.assertEqual(svg.count('class="speedup-label faster"'), 12)
        self.assertEqual(svg.count('class="timing-label"'), 12)
        self.assertEqual(svg.count('data-center-y="197.00"'), 2)
        self.assertEqual(
            len(re.findall(r'class="python-bar"[^>]*height="30"', svg)),
            12,
        )
        self.assertEqual(
            len(re.findall(r'class="julia-bar"[^>]*height="16"', svg)),
            12,
        )
        self.assertRegex(
            svg,
            r'class="workflow-label"[^>]*font-size="15"',
        )
        self.assertRegex(
            svg,
            r'class="speedup-label faster"[^>]*font-size="15"',
        )
        self.assertRegex(
            svg,
            r'class="timing-label"[^>]*font-size="13"',
        )
        self.assertIn("GitHub Actions run 42", svg)
        self.assertIn("0123456789ab", svg)
        self.assertIn("Workflow 12", svg)
        self.assertIn("240.000", svg)
        self.assertIn("120.000", svg)
        self.assertIn("Julia 2.00× faster", svg)
        self.assertIn("P 240.0 · J 120.0 ms", svg)

    def test_speedup_annotations_use_human_readable_direction(self) -> None:
        self.assertEqual(
            _speedup_annotation(2.0, 1.0),
            ("Julia 2.00× faster", "faster"),
        )
        self.assertEqual(
            _speedup_annotation(1.0, 2.0),
            ("Julia 2.00× slower", "slower"),
        )
        self.assertEqual(
            _speedup_annotation(1.0, 1.005),
            ("Julia ≈ parity (1.00×)", "parity"),
        )

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
