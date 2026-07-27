from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path


CI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CI_DIR))

from check_ed_application_catalog import load_catalog, validate_catalog


class EDApplicationCatalogTests(unittest.TestCase):
    def test_committed_catalog_has_fifty_balanced_applications(self) -> None:
        counts = validate_catalog(load_catalog())
        self.assertEqual(len(counts["categories"]), 10)
        self.assertEqual(set(counts["categories"].values()), {5})
        self.assertEqual(sum(counts["coverage"].values()), 50)
        self.assertGreaterEqual(
            counts["automatic_differentiation"]["essential"],
            10,
        )

    def test_missing_gradient_target_is_rejected_for_required_ad(self) -> None:
        catalog = copy.deepcopy(load_catalog())
        application = next(
            item
            for item in catalog["applications"]
            if item["ad_requirement"] == "essential"
        )
        application["gradient_target"] = None
        with self.assertRaisesRegex(ValueError, "gradient_target"):
            validate_catalog(catalog)

    def test_catalog_does_not_claim_new_executable_benchmarks(self) -> None:
        catalog = load_catalog()
        self.assertNotIn("executable", catalog["coverage_definitions"])
        self.assertEqual(
            [item["id"] for item in catalog["applications"]],
            [f"ED-{index:03d}" for index in range(1, 51)],
        )


if __name__ == "__main__":
    unittest.main()
