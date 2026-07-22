# Rust verification preparation

This directory is the test-side foundation for a future Rust QuSpin package.
It is deliberately split into three concerns:

1. `src/api.rs` defines the narrow adapter contract expected from a candidate:
   bases, operator terms, Hamiltonian assembly, a shared `LinearOperator`,
   eigensolvers, Krylov evolution, spectra, subspace fidelity, and Lindblad
   evolution.
2. `src/catalog.rs` ports the exact twelve medium-size Python/Julia paper
   workflows, including their parameters, required capabilities, and physical
   invariants.
3. `src/benchmark.rs` ports the warm-up/sample protocol and emits the same CSV
   fields used by the current timing renderer, with `language=rust`.

`ci/render_timing_report.py` accepts the resulting Rust CSV either alongside
the existing Python/Julia files or as a Python/Rust pair.

The crate does **not** contain a hidden Rust ED implementation. Its current
tests validate the API/harness design, the twelve-case denominator, invariant
enforcement, and CSV compatibility. Once the candidate repository exists, a
private adapter will implement `QuSpinApi` and `WorkflowBackend`; those same
generic tests will then exercise the real crate.

[`full_api_contract.json`](full_api_contract.json) also maps the frozen
64-object / 282-method / 180-attribute migration denominator into the proposed
Rust module boundaries. `ci/check_rust_api_plan.py` prevents either denominator
or namespace coverage from drifting silently. `contract_ready` means only that
the adapter/test boundary exists; it does not mean the Rust implementation or
Python-oracle parity exists.

Run locally with:

```bash
cargo fmt --manifest-path rust/Cargo.toml --check
cargo clippy --manifest-path rust/Cargo.toml --all-targets -- -D warnings
cargo test --manifest-path rust/Cargo.toml
```

The proposed public package surface and the migration map from Julia are in
[`docs/RUST_API_CONTRACT.md`](../docs/RUST_API_CONTRACT.md).
