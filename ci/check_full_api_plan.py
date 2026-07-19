#!/usr/bin/env python3
"""Fail CI if the frozen QuSpin migration denominator becomes inconsistent."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "test" / "full_api_migration_plan.json"


def main() -> None:
    payload = json.loads(PLAN.read_text(encoding="utf-8"))
    counts = payload["counts"]

    assert payload["source_package"] == "quspin"
    assert payload["source_version"] == "1.0.1"
    assert payload["source_snapshot"] == "python-public-api:sha256:47f7f10a6ffcc19a"
    assert counts["objects"] == 64
    assert counts["classes"] == 20
    assert counts["functions"] == 40
    assert counts["values"] == 4
    assert counts["methods_excluding_init"] == 282
    assert counts["attributes"] == 180

    objects = payload["objects"]
    assert len(objects) == counts["objects"]
    assert len({entry["source"] for entry in objects}) == len(objects)
    assert (
        counts["implemented_objects"]
        + counts["partial_objects"]
        + counts["planned_objects"]
        == counts["objects"]
    )

    allowed_statuses = {"planned", "partial", "implemented"}
    for entry in objects:
        assert entry["status"] in allowed_statuses
        tiers = set(entry["required_test_tiers"])
        assert {"surface", "unit", "property"} <= tiers
        if entry["kind"] != "value":
            assert "oracle" in tiers
        if entry["usage"].get("by_kind", {}).get("integration", 0):
            assert "integration" in tiers
        sources = [member["source"] for member in entry["members"]]
        assert len(sources) == len(set(sources))

    print(
        "full API plan:",
        counts["objects"],
        "objects,",
        counts["methods_excluding_init"],
        "methods,",
        counts["attributes"],
        "attributes",
    )


if __name__ == "__main__":
    main()
