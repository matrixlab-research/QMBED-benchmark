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
4. `src/candidate.rs` pins and adapts one public QuSpin.rs commit.
   `src/candidate_workflows.rs` independently constructs every frozen paper
   workflow and returns its named physical observations.

`ci/render_timing_report.py` accepts the resulting Rust CSV either alongside
the existing Python/Julia files or as a Python/Rust pair.

The crate does **not** contain a hidden Rust ED implementation. The private
`WorkflowBackend` calls the pinned public crate and executes all twelve
medium-size workflows. Release-mode CI validates every returned physical
invariant before a benchmark row can be emitted.

[`full_api_contract.json`](full_api_contract.json) also maps the frozen
64-object / 282-method / 180-attribute migration denominator into the proposed
Rust module boundaries. `ci/check_rust_api_plan.py` prevents the denominator,
namespace coverage, or complete-version object mapping from drifting silently.
The low-level adapter and all
twelve workflows run against the pinned candidate, and the first one-warm-up,
five-sample local record is stored under `rust/benchmarks/`. Cross-language
oracle comparison and same-runner hosted timings remain separate gates.

The rewrite task itself is frozen from [`spec_source.json`](spec_source.json)
into the content-addressed bundle under [`specs/quspin`](specs/quspin). Minos
renders the public [`MOTIVATION.md`](taskdoc/MOTIVATION.md),
[`CONTRACT.md`](taskdoc/CONTRACT.md), and [`TESTS.md`](taskdoc/TESTS.md) from
that one bundle. The source and bundle remain private because they contain
held-out anchors; the three generated task documents contain only derivations,
the required Rust surface, and small visible examples.

Those frozen documents specify the accepted 23-symbol paper-workflow core.
They are intentionally not edited as the package grows. The design inputs for
the complete Python-1.0.1 capability denominator are maintained separately in
[`full-taskdoc/MOTIVATION.md`](full-taskdoc/MOTIVATION.md),
[`full-taskdoc/CONTRACT.md`](full-taskdoc/CONTRACT.md), and
[`full-taskdoc/TESTS.md`](full-taskdoc/TESTS.md). The full documents distinguish
surface, semantic, workflow, and scale completion and map every one of the 64
frozen Python objects to a Rust-native boundary.

Run locally with:

```bash
cargo fmt --manifest-path rust/Cargo.toml --check
cargo clippy --manifest-path rust/Cargo.toml --all-targets -- -D warnings
cargo test --manifest-path rust/Cargo.toml
cargo test --release --manifest-path rust/Cargo.toml --test candidate_workflows -- --ignored
QUSPIN_RUST_SAMPLES=5 cargo run --release --manifest-path rust/Cargo.toml --bin paper_benchmark
python ci/freeze_rust_spec.py --check
```

The proposed public package surface and the migration map from Julia are in
[`docs/RUST_API_CONTRACT.md`](../docs/RUST_API_CONTRACT.md).
