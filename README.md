# QuSpinVerify — private verification package

Native Julia package that verifies `matrixlab-research/QuSpin.jl`. Scaffolded
by [minos](https://github.com/kunyuan/minos) `build`; kept **private in the
Matrix Lab organization** so its tests and CI logs are never visible to the
implementing agent.

For the Python-to-Julia QuSpin campaign, the `private-verification` GitHub
Actions workflow first regenerates the frozen observations with a pinned
Python QuSpin commit. It then installs the requested Julia candidate ref and
runs this private Test.jl suite. Python is not present in the candidate job.

```
Project.toml         deps: Test (+ the package under test, added at CI time)
test/runtests.jl     the private suite — ordinary @testset / @test; the answer
                     key. It delegates to namespace/object files under
                     test/full_api/, which mirror the 64-object migration plan.
test/full_api_migration_plan.json
                     frozen 64-object / 282-method / 180-attribute denominator.
ci/runcandidate.jl   develop + test the package the MR proposes.
.github/workflows/verify.yml
                     first reproduces the pinned Python oracle, then runs the
                     private Julia suite against a selected candidate ref.
```

## How a candidate is judged

1. A maintainer dispatches `private-verification` with the public candidate
   repository and an immutable commit SHA or selected ref.
2. The oracle job installs the pinned Python QuSpin commit, regenerates its
   deterministic observations, and compares them structurally against the
   frozen reference. Keys, shapes, integer values, and strings remain exact;
   floating-point values use a `5e-12` relative/absolute tolerance for
   platform-specific BLAS tails.
3. Only after that passes, the Julia job installs the candidate and runs
   `Pkg.test()` with the held-out suite.

Failure logs — including any expected values Test.jl prints — stay in this
private repository. Repository privacy replaces log redaction; no blind
evaluation harness is involved.

## Timing and performance baselines

Every candidate run produces two private timing artifacts:

1. `verification-timing-*` records every API/integration test file's wall,
   compile, recompile, GC, allocation, and lock-conflict measurements. These
   are cold, ordered verification-suite timings and answer which test files
   dominate CI time.
2. `quspin-performance-*` compares the pinned Python QuSpin commit with the
   Julia candidate on the same GitHub Actions runner. It covers basis and
   Hamiltonian construction, dense and CSC matrix-vector action, dense and CSC
   Lanczos, CSC ARPACK partial eigenspectra, dense full diagonalization,
   CSR/DIA kernels, physical momentum-sector construction, matrix-free
   operator construction and action, CSC Krylov time evolution, entanglement
   entropy, two-map 2D general-basis construction, higher-spin sparse
   construction, batched time evolution and entropy, spinful sparse
   construction and matrix-free action, sparse ExpOp vector/grid action,
   diagonal-ensemble fluctuations, and fresh-process package load.
   Dense/CSC/CSR/DIA controlled comparisons are reported separately from
   current-backend comparisons.

Operation benchmarks run three warm-ups followed by 15 samples. Each sample is
adaptively batched to target at least 80 ms, and the raw CSV reports minimum,
p05, p25, median, mean, p75, p95, maximum, standard deviation, and iterations
per sample. Fresh-process package loading uses one warm-up and nine new
processes. Julia allocation medians are recorded separately. BLAS, OpenMP, and
Julia are fixed to one thread. The report explicitly lists Python and Julia
storage formats for every row. A preflight runs outside the timed region;
new-API rows assert dimensions, return shapes, matrix fingerprints, or physical
invariants, while legacy rows retain an exception-free smoke check. Native
Julia CSC construction is measured directly against Python CSC construction.
The suite separately asserts the
stored matrix types, sparse ARPACK residuals, and reduced workflows from
published MBL, quantum-quench, and Floquet spin-chain studies.

The completion gate also keeps independent cases for boson, spinless-fermion,
spinful-fermion, and user-defined Hamiltonians; translation/parity/inversion
sector reconstruction; actual CSR/DIA storage; Hermiticity, particle-number,
and symmetry rejection; matrix-free ARPACK residuals; and sparse Krylov
evolution/Floquet unitarity. These cases intentionally use different lattice
sizes and coefficients from the public package tests. Their numerical anchors
come from the pinned Python commit recorded in `oracle/reference.json`.

The performance job is observational: it reports `Python median / Julia
median` but does not fail a candidate on noisy hosted-runner timing. The CSV
artifacts are retained for 90 days. A regression threshold should only be
introduced after enough runs establish runner variance.

### Workflow-derived ED coverage

The verification suite also contains 21 end-to-end ED scenarios derived from
13 LKM exact-diagonalization workflow families. See
[`docs/ED_WORKFLOW_COVERAGE.md`](docs/ED_WORKFLOW_COVERAGE.md) for the model,
observable, integration size, representative paper scale, and explicit
direct/assisted/proxy boundary of every scenario.
The first paired local paper-scale result is recorded in
[`docs/PAPER_PERFORMANCE_BASELINE.md`](docs/PAPER_PERFORMANCE_BASELINE.md).

Small Julia workflow timings run on every candidate performance job after a
correctness preflight. Select `benchmark_tier=paper` in the manual workflow to
add six paired, paper-shaped Python QuSpin vs Julia measurements. Paper
workflows use lazy per-case construction, one warm-up, five samples, physical
residual/norm/unitarity checks, and preserve raw samples in the uploaded CSV.
They are observational performance evidence rather than a noisy hosted-runner
pass/fail gate.

### Test provenance

The verification inputs have three distinct sources:

1. The API denominator is the frozen
   [64-object migration plan](test/full_api_migration_plan.json) extracted from
   the Python QuSpin surface used for this campaign.
2. Numerical oracle values are regenerated from the pinned upstream
   [QuSpin commit `5bf9e5b`](https://github.com/QuSpin/QuSpin/commit/5bf9e5b266e6d8b70e5cf5973c7c7d59d62e412f)
   and compared with [`oracle/reference.json`](oracle/reference.json).
3. Scientific scenarios start from the 13 exact-diagonalization families
   returned by the [LKM Workflow Families registry](https://lkm.test.bohrium.com/web/en/workflows)
   on 2026-07-20. The resulting 21 deterministic integration cases and their
   direct/assisted/proxy boundaries are listed in
   [`docs/ED_WORKFLOW_COVERAGE.md`](docs/ED_WORKFLOW_COVERAGE.md).

The six timed `paper`-tier workflows are larger paired Python/Julia instances
of that workflow catalog:

| Timed workflow | Workflow/paper source | Scientific kernel retained |
|---|---|---|
| MBL mid-spectrum shift-invert | LKM family `19021`; [Pal and Huse, *The many-body localization transition*](https://arxiv.org/abs/1003.2613) | random-field spin chain, fixed-magnetization sector, interior eigenpairs |
| XXZ Lanczos quench | LKM family `3831`; [Weinberg and Bukov, *QuSpin part I*](https://arxiv.org/abs/1610.03042) | XXZ Hamiltonian and Krylov/Lanczos real-time evolution |
| Floquet full-unitary heating | LKM family `1069`; [Weinberg and Bukov, *QuSpin part I*](https://arxiv.org/abs/1610.03042) | two-step driven spin chain and full Floquet unitary |
| Spinful Hubbard spectrum | LKM family `3127` | finite correlated-fermion cluster and low-energy sparse spectrum |
| Interacting SSH spectrum | LKM family `19000` | half-filled interacting SSH chain and low-energy sparse spectrum |
| Translation-sector XXZ spectrum | LKM family `3831`; [Weinberg and Bukov, *QuSpin part I*](https://arxiv.org/abs/1610.03042) | momentum-sector construction and sparse eigenspectrum |

These are **paper-shaped workflow benchmarks**, not complete paper
reproductions and not timings copied from the papers. Lattice sizes,
coefficients, solver tolerances, and correctness preflights are deterministic
benchmark choices defined in
[`python_paper_workflow_benchmarks.py`](benchmark/python_paper_workflow_benchmarks.py)
and
[`julia_paper_workflow_benchmarks.jl`](benchmark/julia_paper_workflow_benchmarks.jl).

### Latest six-workflow timing chart

The grouped bar chart below is updated automatically after a successful
`paper`-tier run dispatched from `main`. It reports the latest same-runner
median wall time for Python QuSpin and the selected Julia candidate; lower is
better. PR and non-`main` runs still upload the SVG as an Actions artifact but
do not rewrite this README snapshot.

![Python and Julia median timings for six paper workflows](docs/paper-workflow-timings.svg)

## Extending the suite

Each namespace has a directory under `test/full_api/`; add one file per public
object or tightly coupled helper family, then include it from `test/runtests.jl`.
The held-out goldens are a floor, not a ceiling: add randomized property checks
and real multi-call integration tests. `ci/check_full_api_plan.py` prevents the
source denominator from silently shrinking.
