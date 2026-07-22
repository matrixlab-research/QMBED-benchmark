#!/usr/bin/env python3
"""Freeze the authored Rust rewrite spec with Minos-compatible hashes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
from typing import Any


def canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def content_hash(value: Any) -> str:
    digest = hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


def freeze(source: Path, output_root: Path) -> Path:
    bundle = json.loads(source.read_text(encoding="utf-8"))
    if bundle.get("target_language") != "rust":
        raise SystemExit("Rust source bundle must set target_language=rust")
    if not bundle.get("specs"):
        raise SystemExit("Rust source bundle has no specs")

    symbols: set[str] = set()
    for spec in bundle["specs"]:
        if spec.get("package") != bundle["package"]:
            raise SystemExit(f"package mismatch in spec {spec.get('symbol', '?')}")
        symbol = spec.get("symbol", "")
        if not symbol or symbol in symbols:
            raise SystemExit(f"missing or duplicate Rust spec symbol: {symbol!r}")
        symbols.add(symbol)
        motivation = spec.get("motivation") or {}
        contract = spec.get("contract") or {}
        anchors = spec.get("anchors") or []
        properties = spec.get("properties") or []
        visible = sum(not anchor.get("held_out", False) for anchor in anchors)
        held_out = len(anchors) - visible
        pinning = len(properties) + 2 * held_out
        if not (motivation.get("question") or motivation.get("definition")):
            raise SystemExit(f"thin Rust spec {symbol}: missing motivation")
        if not (
            contract.get("comparison")
            or contract.get("atol") is not None
            or contract.get("rtol") is not None
        ):
            raise SystemExit(f"thin Rust spec {symbol}: missing comparison")
        if visible < 1 or held_out < 1 or pinning < 10:
            raise SystemExit(
                f"thin Rust spec {symbol}: visible={visible}, held_out={held_out}, "
                f"properties={len(properties)}, pinning={pinning}"
            )
        spec.pop("spec_hash", None)
        spec["spec_hash"] = content_hash(spec)

    payload = {
        "package": bundle["package"],
        "version": bundle["version"],
        "spec_hashes": sorted(spec["spec_hash"] for spec in bundle["specs"]),
    }
    bundle["bundle_hash"] = content_hash(payload)

    package_dir = output_root / bundle["package"]
    package_dir.mkdir(parents=True, exist_ok=True)
    short_hash = bundle["bundle_hash"].split(":", maxsplit=1)[1][:12]
    artifact_name = f"{bundle['version']}.{short_hash}.json"
    artifact = package_dir / artifact_name
    artifact.write_text(
        json.dumps(bundle, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    latest = {
        "package": bundle["package"],
        "version": bundle["version"],
        "target_language": "rust",
        "bundle_hash": bundle["bundle_hash"],
        "artifact": artifact_name,
        "snapshot": bundle.get("snapshot", ""),
        "frozen_at": bundle.get("frozen_at", ""),
    }
    (package_dir / "latest.json").write_text(
        json.dumps(latest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return artifact


def check(source: Path, committed_root: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="quspin-rust-spec-") as directory:
        generated_root = Path(directory)
        generated_artifact = freeze(source, generated_root)
        generated_latest = generated_artifact.parent / "latest.json"

        committed_package = committed_root / "quspin"
        committed_latest = committed_package / "latest.json"
        if not committed_latest.is_file():
            raise SystemExit(f"missing committed pointer: {committed_latest}")
        latest = json.loads(committed_latest.read_text(encoding="utf-8"))
        committed_artifact = committed_package / latest["artifact"]
        if not committed_artifact.is_file():
            raise SystemExit(f"missing committed artifact: {committed_artifact}")

        if generated_latest.read_bytes() != committed_latest.read_bytes():
            raise SystemExit("rust/specs/quspin/latest.json is stale; run freeze_rust_spec.py")
        if generated_artifact.read_bytes() != committed_artifact.read_bytes():
            raise SystemExit(f"{committed_artifact} is stale; run freeze_rust_spec.py")
        extras = sorted(
            path.name
            for path in committed_package.glob("*.json")
            if path.name not in {"latest.json", committed_artifact.name}
        )
        if extras:
            raise SystemExit(f"stale frozen Rust spec artifacts remain: {', '.join(extras)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path("rust/spec_source.json"))
    parser.add_argument("--out", type=Path, default=Path("rust/specs"))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        check(args.source, args.out)
        print("Rust SpecBundle is current")
    else:
        artifact = freeze(args.source, args.out)
        print(artifact)


if __name__ == "__main__":
    main()
