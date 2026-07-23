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

| Workflow | First baseline (s) | Complete candidate (s) | Universal backends (s) | Candidate → universal |
|---|---:|---:|---:|---:|
| MBL shift-invert | 0.420 | 0.427 | 0.403 | 1.06× |
| XXZ quench | 0.056 | 0.062 | 0.035 | 1.79× |
| Floquet full unitary | 1.094 | 1.102 | 1.080 | 1.02× |
| Spinful Hubbard spectrum | 0.113 | 0.114 | 0.070 | 1.63× |
| Interacting SSH spectrum | 0.295 | 0.303 | 0.187 | 1.62× |
| Translation-sector XXZ | 0.038 | 0.045 | 0.028 | 1.57× |
| TFIM fidelity scan | 5.426 | 5.471 | 2.366 | 2.31× |
| PXP revival | 0.672 | 0.689 | 0.313 | 2.20× |
| Bose-Hubbard quench | 0.171 | 0.170 | 0.083 | 2.05× |
| Spinful-Hubbard current | 2.719 | 2.721 | 1.644 | 1.66× |
| CoNb2O6 response | 4.259 | 4.360 | 2.058 | 2.12× |
| Particle-addition spectrum | 1.274 | 1.321 | 0.725 | 1.82× |

The new record used `rustc 1.97.1`, one validated warm-up, and five measured
samples on the same local macOS arm64 host. Relative to the complete-candidate
record, the geometric-mean speedup is 1.69× and the sum of the twelve medians
falls from 16.79 s to 8.99 s, a 1.87× end-to-end improvement. The complete
candidate used `rustc 1.85.0`, so this remains an observational implementation
comparison rather than a compiler-controlled microbenchmark; the original
`rustc 1.97.1` record is retained as an additional continuity anchor.

Regenerate the CSV with:

```bash
QUSPIN_RUST_WARMUPS=1 QUSPIN_RUST_SAMPLES=5 \
  cargo run --release --manifest-path rust/Cargo.toml --bin paper_benchmark
```
