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

| Workflow | First baseline (s) | Complete candidate (s) |
|---|---:|---:|
| MBL shift-invert | 0.420 | 0.427 |
| XXZ quench | 0.056 | 0.062 |
| Floquet full unitary | 1.094 | 1.102 |
| Spinful Hubbard spectrum | 0.113 | 0.114 |
| Interacting SSH spectrum | 0.295 | 0.303 |
| Translation-sector XXZ | 0.038 | 0.045 |
| TFIM fidelity scan | 5.426 | 5.471 |
| PXP revival | 0.672 | 0.689 |
| Bose-Hubbard quench | 0.171 | 0.170 |
| Spinful-Hubbard current | 2.719 | 2.721 |
| CoNb2O6 response | 4.259 | 4.360 |
| Particle-addition spectrum | 1.274 | 1.321 |

Because the two stored records use different compiler versions, the table is
an observational continuity check rather than a controlled compiler A/B. Ten
of twelve medians remain within 4% of the first baseline; the two
transition-dominated XXZ cases differ by 11% and 19%, respectively.

Regenerate the CSV with:

```bash
QUSPIN_RUST_WARMUPS=1 QUSPIN_RUST_SAMPLES=5 \
  cargo run --release --manifest-path rust/Cargo.toml --bin paper_benchmark
```
