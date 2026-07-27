# Native AD acceptance benchmarks

The AD suite validates a scientific derivative at the end-to-end workflow
boundary. It does not time an isolated dual-number operation and does not
differentiate through sparse assembly or Lanczos instructions.

For every workflow, QMBED computes the isolated ground-state energy gradient
with the Hellmann--Feynman rule. The independent oracle evaluates a centered
finite difference using the same operator family, eigensolver target,
tolerance, and reusable eigensolver workspace.

## Workflow coverage

| Case | Scientific class | Derived from | Catalog / LKM family |
|---|---|---|---|
| `ad_tfim_ground_energy` | transverse-field Ising model | TFIM fidelity workflow | ED-003, ED-041 / 30690, 23806, 3831 |
| `ad_xxz_sector_ground_energy` | fixed-magnetization XXZ chain | XXZ quench workflow | ED-001 / 3831, 4905 |
| `ad_j1_j2_frustrated_ground_energy` | frustrated J1-J2 chain | frustrated-spin workflow | ED-002 / 23806, 3831 |
| `ad_translation_xxz_ground_energy` | translation-reduced XXZ chain | translation-sector workflow | ED-001 / 3831, 4905 |
| `ad_interacting_ssh_ground_energy` | interacting SSH chain | interacting SSH workflow | ED-016 / 19000, 16282 |
| `ad_disordered_tv_ground_energy` | disordered spinless t-V model | MBL workflow | ED-034 / 19021, 3831 |
| `ad_triangular_fermion_ground_energy` | triangular-lattice fermions | particle-addition proxy | ED-018 / 26341 |
| `ad_hubbard_ground_energy` | spinful Hubbard model | Hubbard spectrum workflow | ED-006 / 32600 |
| `ad_ionic_hubbard_ground_energy` | ionic Hubbard model | Hubbard transport/phase workflow | ED-005, ED-006 / 8704, 3127, 32600 |
| `ad_bose_hubbard_ground_energy` | Bose-Hubbard model | Mott-quench workflow | paper workflow; not separately catalogued |
| `ad_trapped_bose_hubbard_ground_energy` | trapped Bose-Hubbard model | inhomogeneous boson workflow | paper workflow; not separately catalogued |
| `ad_pxp_detuning_ground_energy` | constrained PXP basis | PXP revival workflow | ED-009 / 22708 |

This is six broad capability classes: unrestricted and symmetry-reduced spins,
spinless fermions, spinful fermions, bosons, non-chain geometry, and a
callback-defined constrained basis.

ED and LKM identifiers refer to the reviewed
[50-application catalog](LKM_ED_APPLICATION_CATALOG.md). Paper citations and
the precise direct/proxy relationship of the source workflows are listed in
the repository [README](../README.md#test-provenance). The AD cases reuse their
scientific model families but define a new differentiable target; they do not
claim to reproduce the cited papers.

## Acceptance gates

Each CSV row records:

- maximum absolute and componentwise relative gradient error;
- ground-state residual and spectral gap;
- analytic and finite-difference wall time;
- parameter count and eigensolve counts;
- finite-difference time divided by analytic time.

A case passes when either the absolute error is at most `2e-5` or the relative
error is at most `2e-4`, the residual is at most `1e-8`, and the ground-state
gap is resolved. Analytic AD must use one eigensolve; central differences use
two eigensolves per parameter. The measured geometric-mean speedup across all
twelve cases must exceed one.

The thresholds are committed here before CI produces the measured result. They
are not adjusted per workflow.

## CI profiles

Pull requests run all twelve cases with smaller, nontrivial Hilbert spaces and
one timing sample. A scheduled Monday job and manual `paper` dispatch use
medium ED sizes, one warm-up, and three samples. Both upload
`ad-benchmark.csv`; timing is the median and the execution order alternates to
reduce first/second-run bias.

Run locally:

```bash
QMBED_AD_SCALE=ci QMBED_AD_WARMUPS=0 QMBED_AD_SAMPLES=1 \
  cargo run --release --manifest-path rust/Cargo.toml --bin ad_benchmark \
  > ad-benchmark.csv
python ci/check_ad_benchmark.py ad-benchmark.csv
```

The committed 2026-07-28 paper-size record is
[`rust/benchmarks/native_ad_paper_2026-07-28.csv`](../rust/benchmarks/native_ad_paper_2026-07-28.csv).
With one warm-up and three samples, its geometric-mean speedup is 7.44×
(3.92×–13.88× across individual cases) for dimensions 495–16,384. The
CI-size five-sample record is retained separately so the fast gate and
paper-size evidence are not conflated.

## Current boundary

These workflows validate operator JVP/VJP and isolated ground-state energy
gradients. They do not imply that eigenvector, degenerate-subspace, time
evolution, thermal-trace, Floquet, or response-function rules are implemented.
Those remain separate native-AD capability work.
