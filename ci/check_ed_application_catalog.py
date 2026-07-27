from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = ROOT / "benchmark" / "catalog" / "lkm_ed_applications.json"

EXPECTED_CATEGORIES = {
    "equilibrium_and_phase_diagrams",
    "dynamics_and_transport",
    "response_and_spectroscopy",
    "topology_and_entanglement",
    "open_systems_and_quantum_optics",
    "cross_domain_ed",
    "thermodynamics_localization_and_certification",
    "metrology_and_control",
    "differentiable_inverse_problems",
    "differentiable_numerical_kernels",
}
VALID_COVERAGE = {"direct", "partial", "beyond"}
VALID_AD_REQUIREMENT = {"none", "helpful", "essential"}
REQUIRED_APPLICATION_FIELDS = {
    "id",
    "category",
    "title",
    "lkm_family_ids",
    "qmbed_coverage",
    "ad_requirement",
    "gradient_target",
    "benchmark_shape",
    "primary_observable",
    "validation",
    "capability_gaps",
    "evidence_dois",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_catalog(path: Path = CATALOG_PATH) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def validate_catalog(catalog: dict[str, Any]) -> dict[str, Counter[str]]:
    _require(catalog.get("schema_version") == 1, "schema_version must be 1")
    source = catalog.get("source")
    _require(isinstance(source, dict), "source must be an object")
    _require(
        source.get("workflow_index") == "https://lkm.bohrium.com/web/zh/workflows",
        "workflow_index must point to the LKM Workflow Families page",
    )
    _require(
        source.get("workflow_queries", 0) >= 90,
        "the search manifest must retain the systematic query count",
    )
    _require(
        source.get("unique_families_screened", 0) >= 500,
        "the catalog must record the screened-family denominator",
    )

    applications = catalog.get("applications")
    _require(isinstance(applications, list), "applications must be an array")
    _require(len(applications) == 50, "catalog must contain exactly 50 applications")

    expected_ids = [f"ED-{index:03d}" for index in range(1, 51)]
    actual_ids = [application.get("id") for application in applications]
    _require(actual_ids == expected_ids, "application IDs must be ED-001 through ED-050")

    for application in applications:
        app_id = application["id"]
        missing = REQUIRED_APPLICATION_FIELDS - set(application)
        _require(not missing, f"{app_id} is missing fields: {sorted(missing)}")
        _require(
            application["category"] in EXPECTED_CATEGORIES,
            f"{app_id} has an unknown category",
        )
        _require(
            application["qmbed_coverage"] in VALID_COVERAGE,
            f"{app_id} has an invalid qmbed_coverage",
        )
        _require(
            application["ad_requirement"] in VALID_AD_REQUIREMENT,
            f"{app_id} has an invalid ad_requirement",
        )
        _require(
            isinstance(application["lkm_family_ids"], list)
            and application["lkm_family_ids"]
            and all(
                isinstance(cluster_id, int) and cluster_id > 0
                for cluster_id in application["lkm_family_ids"]
            ),
            f"{app_id} must cite at least one positive LKM family ID",
        )
        for field in ("title", "benchmark_shape", "primary_observable", "validation"):
            _require(
                isinstance(application[field], str) and application[field].strip(),
                f"{app_id}.{field} must be non-empty",
            )
        _require(
            isinstance(application["capability_gaps"], list)
            and all(
                isinstance(gap, str) and gap.strip()
                for gap in application["capability_gaps"]
            ),
            f"{app_id}.capability_gaps must be an array of strings",
        )
        _require(
            isinstance(application["evidence_dois"], list)
            and all(
                isinstance(doi, str) and doi.startswith("10.")
                for doi in application["evidence_dois"]
            ),
            f"{app_id}.evidence_dois must contain DOI strings",
        )
        if application["ad_requirement"] == "essential":
            _require(
                isinstance(application["gradient_target"], str)
                and application["gradient_target"].strip(),
                f"{app_id} requires a concrete gradient_target",
            )
            _require(
                application["qmbed_coverage"] == "beyond",
                f"{app_id} cannot be marked covered while required AD is absent",
            )
            _require(
                application["evidence_dois"],
                f"{app_id} requires at least one method or application DOI",
            )

    category_counts = Counter(
        application["category"] for application in applications
    )
    _require(
        set(category_counts) == EXPECTED_CATEGORIES,
        "catalog must cover all ten categories",
    )
    _require(
        set(category_counts.values()) == {5},
        "each category must contain exactly five applications",
    )

    coverage_counts = Counter(
        application["qmbed_coverage"] for application in applications
    )
    ad_counts = Counter(
        application["ad_requirement"] for application in applications
    )
    _require(
        coverage_counts["beyond"] >= 20,
        "the frontier catalog must retain at least 20 beyond-QMBED applications",
    )
    _require(
        ad_counts["essential"] >= 10,
        "the catalog must retain at least 10 AD-essential applications",
    )
    return {
        "categories": category_counts,
        "coverage": coverage_counts,
        "automatic_differentiation": ad_counts,
    }


def main() -> None:
    counts = validate_catalog(load_catalog())
    print(
        "validated 50 ED applications: "
        f"{len(counts['categories'])} categories, "
        f"{counts['coverage']['direct']} direct, "
        f"{counts['coverage']['partial']} partial, "
        f"{counts['coverage']['beyond']} beyond QMBED, "
        f"{counts['automatic_differentiation']['essential']} AD-essential"
    )


if __name__ == "__main__":
    main()
