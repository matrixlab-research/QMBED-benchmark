#!/usr/bin/env python3
"""Keep the Rust contract tied to the frozen 64-object API denominator."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PLAN = ROOT / "test" / "full_api_migration_plan.json"
RUST_CONTRACT = ROOT / "rust" / "full_api_contract.json"


def validate(source_path: Path, contract_path: Path, root: Path) -> str:
    source = json.loads(source_path.read_text())
    contract = json.loads(contract_path.read_text())

    expected_counts = contract["counts"]
    for key in ("objects", "methods_excluding_init", "attributes"):
        actual = source["counts"][key]
        expected = expected_counts[key]
        if actual != expected:
            raise ValueError(
                f"Rust denominator drift for {key}: expected {expected}, got {actual}"
            )

    source_namespaces = Counter(item["namespace"] for item in source["objects"])
    mapped_namespaces = contract["namespaces"]
    if set(source_namespaces) != set(mapped_namespaces):
        raise ValueError(
            "Rust namespace map differs from source denominator: "
            f"source={sorted(source_namespaces)}, "
            f"mapped={sorted(mapped_namespaces)}"
        )

    allowed_states = set(contract["activation_states"])
    for namespace, actual_count in source_namespaces.items():
        mapping = mapped_namespaces[namespace]
        if mapping["objects"] != actual_count:
            raise ValueError(
                f"Rust namespace {namespace} maps {mapping['objects']} objects; "
                f"source has {actual_count}"
            )
        if not mapping["rust_modules"]:
            raise ValueError(f"Rust namespace {namespace} has no public module")
        if mapping["status"] not in allowed_states:
            raise ValueError(
                f"Rust namespace {namespace} has unknown status {mapping['status']}"
            )
        for relative_path in mapping["verification_contracts"]:
            if not (root / relative_path).is_file():
                raise ValueError(
                    f"Rust namespace {namespace} references missing contract "
                    f"{relative_path}"
                )

    allowed_tiers = {"surface", "unit", "property", "integration", "oracle"}
    for item in source["objects"]:
        declared_tiers = set(item["required_test_tiers"])
        if not declared_tiers:
            raise ValueError(
                f"source object {item['source']} has no required test tier"
            )
        unknown_tiers = declared_tiers.difference(allowed_tiers)
        if unknown_tiers:
            raise ValueError(
                f"source object {item['source']} has unknown tiers "
                f"{sorted(unknown_tiers)}"
            )

    return (
        "Rust API contract: "
        f"{expected_counts['objects']} objects, "
        f"{expected_counts['methods_excluding_init']} methods, "
        f"{expected_counts['attributes']} attributes; "
        + ", ".join(
            f"{namespace}={count}"
            for namespace, count in sorted(source_namespaces.items())
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE_PLAN)
    parser.add_argument("--contract", type=Path, default=RUST_CONTRACT)
    args = parser.parse_args()
    print(validate(args.source, args.contract, ROOT))


if __name__ == "__main__":
    main()
