# Rust verification preparation

This directory is the private test-side foundation for QuSpin.rs. It is
deliberately split into four concerns:

1. `src/api.rs` defines the narrow adapter contract expected from a candidate:
   bases, operator terms, Hamiltonian assembly, a shared `LinearOperator`,
   eigensolvers, Krylov evolution, spectra, subspace fidelity, and Lindblad
   evolution.
2. `src/catalog.rs` ports the exact twelve medium-size Python/Julia paper
   workflows, including their parameters, required capabilities, and physical
   invariants.
3. `src/benchmark.rs` ports the warm-up/sample protocol and emits the same CSV
   fields used by the current timing renderer, with `language=rust`.
4. `src/candidate.rs` pins and adapts one public QuSpin.rs commit. The adapter
   smoke test exercises real basis construction, universal assembly, and
   `eigsh`; it is no longer a contract-only mock.

`ci/render_timing_report.py` accepts the resulting Rust CSV either alongside
the existing Python/Julia files or as a Python/Rust pair.

The crate does **not** contain a hidden Rust ED implementation. Its current
tests validate the API/harness design, the twelve-case denominator, invariant
enforcement, CSV compatibility, and the low-level adapter against the public
crate. `WorkflowBackend` is intentionally still inactive until each paper case
has a real implementation and passes its physical invariants.

[`full_api_contract.json`](full_api_contract.json) also maps the frozen
64-object / 282-method / 180-attribute migration denominator into the proposed
Rust module boundaries. `ci/check_rust_api_plan.py` prevents either denominator
or namespace coverage from drifting silently. `adapter_compiles` means only
that the pinned public crate crosses the private low-level boundary; it does
not mean Python-oracle parity or paper-workflow completion.

The rewrite task itself is frozen from [`spec_source.json`](spec_source.json)
into the content-addressed bundle under [`specs/quspin`](specs/quspin). Minos
renders the public [`MOTIVATION.md`](taskdoc/MOTIVATION.md),
[`CONTRACT.md`](taskdoc/CONTRACT.md), and [`TESTS.md`](taskdoc/TESTS.md) from
that one bundle. The source and bundle remain private because they contain
held-out anchors; the three generated task documents contain only derivations,
the required Rust surface, and small visible examples.

Run locally with:

```bash
cargo fmt --manifest-path rust/Cargo.toml --check
cargo clippy --manifest-path rust/Cargo.toml --all-targets -- -D warnings
cargo test --manifest-path rust/Cargo.toml
python ci/freeze_rust_spec.py --check
```

The proposed public package surface and the migration map from Julia are in
[`docs/RUST_API_CONTRACT.md`](../docs/RUST_API_CONTRACT.md).
