# QuSpinVerify — private verification package

Native Julia package that verifies `kunyuan/QuSpin.jl`. Scaffolded by
[minos](https://github.com/kunyuan/minos) `build`; kept **private on the
maintainer's account** so its tests and CI logs are never visible to the
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
   Lanczos, dense full diagonalization, current-backend time evolution,
   entanglement entropy, and fresh-process package load. Dense/CSC controlled
   comparisons are reported separately from current-backend comparisons.

Operation benchmarks run three warm-ups followed by 15 samples. Each sample is
adaptively batched to target at least 80 ms, and the raw CSV reports minimum,
p05, p25, median, mean, p75, p95, maximum, standard deviation, and iterations
per sample. Fresh-process package loading uses one warm-up and nine new
processes. Julia allocation medians are recorded separately. BLAS, OpenMP, and
Julia are fixed to one thread. The report explicitly lists Python and Julia
storage formats for every row. Julia's current lack of native sparse
Hamiltonian construction is reported as an unsupported capability rather than
timed against Python sparse construction as though the formats were equal.

The performance job is observational: it reports `Python median / Julia
median` but does not fail a candidate on noisy hosted-runner timing. The CSV
artifacts are retained for 90 days. A regression threshold should only be
introduced after enough runs establish runner variance.

## Extending the suite

Each namespace has a directory under `test/full_api/`; add one file per public
object or tightly coupled helper family, then include it from `test/runtests.jl`.
The held-out goldens are a floor, not a ceiling: add randomized property checks
and real multi-call integration tests. `ci/check_full_api_plan.py` prevents the
source denominator from silently shrinking.
