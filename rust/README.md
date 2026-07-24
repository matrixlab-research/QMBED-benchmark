# QMBED Rust benchmark adapter

This directory is the public, independent benchmark-side foundation for QMBED. It is
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
4. `src/candidate.rs` pins and adapts one public QMBED commit.
   `src/candidate_workflows.rs` independently constructs every frozen paper
   workflow and returns its named physical observations.
5. `tests/full_api_oracles.rs` exercises independent full-package semantics that
   are not covered by the original 23-symbol adapter: general symmetry
   projectors, higher-spin and Majorana algebra, wide and branching bases,
   photon sectors, dynamic/parameterized operators, mixed-state measurements,
   general exponential action, and dense/sparse archives.
6. `tests/full_api_scale.rs` checks multiple-size combinatorial enumeration,
   two-map symmetry intersections, direct cross-sector action, sparse algebra
   memory, arbitrary subsystems, and state widths above `u128`.

`ci/render_timing_report.py` accepts the resulting Rust CSV either alongside
the existing Python/Julia files or as a Python/Rust pair.

The crate does **not** contain a second Rust ED implementation. Its
`WorkflowBackend` calls the pinned public crate and executes all twelve
medium-size workflows. Release-mode CI validates every returned physical
invariant before a benchmark row can be emitted.

[`full_api_contract.json`](full_api_contract.json) also maps the frozen
64-object / 282-method / 180-attribute migration denominator into the proposed
Rust module boundaries. `ci/check_rust_api_plan.py` prevents the denominator,
namespace coverage, or complete-version object mapping from drifting silently.
The low-level adapter, independent full-API oracles, structural scale gates, and
all twelve workflows run against the pinned candidate. Current completion
evidence is summarized in [`FULL_API_STATUS.md`](FULL_API_STATUS.md). The first
baseline, complete candidate, universal assembler, and numerical-backend
one-warm-up/five-sample local records are stored under `rust/benchmarks/`.
Same-runner hosted Python/Julia/Rust timing remains a separate observational
gate.

The rewrite task itself is frozen from [`spec_source.json`](spec_source.json)
into the content-addressed bundle under [`specs/quspin`](specs/quspin). Minos
renders [`MOTIVATION.md`](taskdoc/MOTIVATION.md),
[`CONTRACT.md`](taskdoc/CONTRACT.md), and [`TESTS.md`](taskdoc/TESTS.md) from
that one bundle. The source, bundle, and generated documents are now public.
Their `held_out` fields are frozen provenance labels from the original
clean-room campaign, not claims about current repository visibility.

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
cargo test --manifest-path rust/Cargo.toml --test full_api_scale
cargo test --release --manifest-path rust/Cargo.toml --test candidate_workflows -- --ignored
QMBED_RUST_SAMPLES=5 cargo run --release --manifest-path rust/Cargo.toml --bin paper_benchmark
python ci/freeze_rust_spec.py --check
```

The proposed public package surface and the migration map from Julia are in
[`docs/RUST_API_CONTRACT.md`](../docs/RUST_API_CONTRACT.md).
