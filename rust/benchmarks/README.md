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

| Workflow | Rust median (s) |
|---|---:|
| MBL shift-invert | 0.420 |
| XXZ quench | 0.056 |
| Floquet full unitary | 1.094 |
| Spinful Hubbard spectrum | 0.113 |
| Interacting SSH spectrum | 0.295 |
| Translation-sector XXZ | 0.038 |
| TFIM fidelity scan | 5.426 |
| PXP revival | 0.672 |
| Bose-Hubbard quench | 0.171 |
| Spinful-Hubbard current | 2.719 |
| CoNb2O6 response | 4.259 |
| Particle-addition spectrum | 1.274 |

Regenerate the CSV with:

```bash
QUSPIN_RUST_WARMUPS=1 QUSPIN_RUST_SAMPLES=5 \
  cargo run --release --manifest-path rust/Cargo.toml --bin paper_benchmark
```
