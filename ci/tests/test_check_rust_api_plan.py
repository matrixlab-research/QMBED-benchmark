from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


CI_DIR = Path(__file__).resolve().parents[1]
ROOT = CI_DIR.parent
sys.path.insert(0, str(CI_DIR))

from check_rust_api_plan import RUST_CONTRACT, SOURCE_PLAN, validate


class RustApiPlanTests(unittest.TestCase):
    def test_repository_contract_covers_the_frozen_denominator(self) -> None:
        summary = validate(SOURCE_PLAN, RUST_CONTRACT, ROOT)
        self.assertIn("64 objects", summary)
        self.assertIn("282 methods", summary)
        self.assertIn("180 attributes", summary)

    def test_denominator_drift_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            changed = Path(temporary) / "contract.json"
            contract = json.loads(RUST_CONTRACT.read_text())
            contract["counts"]["objects"] = 63
            changed.write_text(json.dumps(contract))
            with self.assertRaisesRegex(ValueError, "denominator drift"):
                validate(SOURCE_PLAN, changed, ROOT)


if __name__ == "__main__":
    unittest.main()
