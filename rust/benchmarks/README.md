# Rust paper-workflow benchmark

`paper_workflows_2026-07-22_local.csv` is the first complete Rust record for
the frozen twelve-workflow denominator. It pins QuSpin.rs commit `93447a4`,
runs one correctness-checked warm-up and five end-to-end samples per case, and
includes basis construction, Hamiltonian assembly, solver/evolution, and the
final observable. Every sample is accepted only after the case's named
physical invariants pass.

The run used `rustc 1.97.1` on a local macOS arm64 machine. These values are a
reproducible local baseline, not a hosted-runner comparison against the
existing Python and Julia records.

`paper_workflows_2026-07-23_local.csv` is the corresponding five-sample record
for the complete candidate at QuSpin.rs commit
`5656fa5482a5bdeb2001ceec82a778573efc717b`, run with Rust 1.85.0 on the same
local macOS arm64 machine. This candidate uses one universal streamed
transition interface: every basis feeds the same assembler, zero to two local
destinations remain inline, and symmetry reduction is fused with target-index
lookup. It does not contain workflow-specific assemblers.

`paper_workflows_2026-07-23_universal_backends.csv` records the merged
structure-aware implementation at QuSpin.rs commit
`42887668e1d6a7373185dc99f252bc080b004be2`. It uses the same public workflow
and universal transition contract, but compiles operator symbols once,
coalesces sparse entries within each source column, selects structural state
indexers, and uses real Lanczos arithmetic when the stored Hamiltonian is
exactly real. None of these paths dispatches on a workflow or model name.

`paper_workflows_2026-07-23_backend_module.csv` pins QuSpin.rs commit
`85a3f5d44b2f5b5bce58e0c35b3e6d31148d43ce`. Dense decomposition, matrix
products, and sparse shifted factorization now cross one internal backend
boundary. Real CSC shift-invert retains real arithmetic through factorization
and Lanczos; the algorithm is selected from operator capabilities, never a
workflow name.

| Workflow | First baseline (s) | Complete candidate (s) | Universal backends (s) | Backend module (s) | Universal → backend |
|---|---:|---:|---:|---:|---:|
| MBL shift-invert | 0.420 | 0.427 | 0.403 | 0.114 | **3.52×** |
| XXZ quench | 0.056 | 0.062 | 0.035 | 0.036 | 0.97× |
| Floquet full unitary | 1.094 | 1.102 | 1.080 | 1.078 | 1.00× |
| Spinful Hubbard spectrum | 0.113 | 0.114 | 0.070 | 0.071 | 0.98× |
| Interacting SSH spectrum | 0.295 | 0.303 | 0.187 | 0.199 | 0.94× |
| Translation-sector XXZ | 0.038 | 0.045 | 0.028 | 0.030 | 0.95× |
| TFIM fidelity scan | 5.426 | 5.471 | 2.366 | 2.407 | 0.98× |
| PXP revival | 0.672 | 0.689 | 0.313 | 0.311 | 1.01× |
| Bose-Hubbard quench | 0.171 | 0.170 | 0.083 | 0.083 | 1.01× |
| Spinful-Hubbard current | 2.719 | 2.721 | 1.644 | 1.661 | 0.99× |
| CoNb2O6 response | 4.259 | 4.360 | 2.058 | 2.098 | 0.98× |
| Particle-addition spectrum | 1.274 | 1.321 | 0.725 | 0.703 | 1.03× |

The new record used `rustc 1.97.1`, one validated warm-up, and five measured
samples on the same local macOS arm64 host. Relative to the complete-candidate
record, the backend-module geometric-mean speedup is 1.85× and the sum of the
twelve medians falls from 16.79 s to 8.79 s, a 1.91× end-to-end improvement.
Relative to the preceding universal-backend record, the geometric mean
improves 1.10× and the sum improves 1.02×; the dominant change is a 3.52× MBL
shift-invert gain. The complete
candidate used `rustc 1.85.0`, so this remains an observational implementation
comparison rather than a compiler-controlled microbenchmark; the original
`rustc 1.97.1` record is retained as an additional continuity anchor.

Regenerate the CSV with:

```bash
QMBED_RUST_WARMUPS=1 QMBED_RUST_SAMPLES=5 \
  cargo run --release --manifest-path rust/Cargo.toml --bin paper_benchmark
```

## Latest compatibility and universal-Krylov A/B

`qmbed_compat_performance_ab_2026-07-26.csv` compares the pre-compatibility
QMBED baseline `58360f5ffed845d0f2fbd3649ad7588420edbc83` with the merged
compatibility candidate `bc7159e6370b3a9c55eede50cf99c8a7df8b9847`. Both sides used the same Apple
M3 Max host, macOS arm64, Rust 1.97.1, a release build, one validated warm-up,
and five complete samples per workflow.

`paper_workflows_2026-07-26_qmbed_latest.csv` is a fresh five-sample run made
after Cargo resolved that exact merged revision from GitHub. It preserves every
raw sample, physical invariant, and evolution-work counter.

| Workflow | Baseline (s) | Candidate (s) | Candidate / baseline |
|---|---:|---:|---:|
| MBL shift-invert | 0.132 | 0.128 | 0.974 |
| XXZ quench | 0.036 | 0.033 | 0.923 |
| Floquet full unitary | 1.126 | 1.098 | 0.975 |
| Spinful Hubbard spectrum | 0.072 | 0.074 | 1.027 |
| Interacting SSH spectrum | 0.201 | 0.200 | 0.993 |
| Translation-sector XXZ | 0.032 | 0.031 | 0.987 |
| TFIM fidelity scan | 2.512 | 2.477 | 0.986 |
| PXP revival | 0.327 | 0.446 | 1.365 |
| Bose-Hubbard quench | 0.084 | 0.319 | 3.814 |
| Spinful-Hubbard current | 1.709 | 1.752 | 1.025 |
| CoNb2O6 response | 2.108 | 1.995 | 0.947 |
| Particle-addition spectrum | 0.735 | 0.719 | 0.979 |

The geometric-mean candidate/baseline ratio is **1.129**, below the 1.15
acceptance ceiling. The summed medians are 9.272 s versus 9.071 s
(`1.022x`). Ten workflows remain between `0.923x` and `1.027x`.

PXP and Bose-Hubbard intentionally perform more work because the baseline used
one 100-vector Krylov projection for the complete time grid and produced
order-one state errors (0.596 and 0.559). The adaptive candidate reaches errors
of 0 and `4.84e-10` against a tighter `1e-12` run. Its observable algorithmic
work is respectively 2 projections / 200 matrix-vector products and
8 projections / 800 matrix-vector products. The first projection in each case
uses the universal real-operator/real-state basis; no workflow or model name is
used for dispatch.
