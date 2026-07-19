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
   deterministic observations, and diffs them against the frozen reference.
3. Only after that passes, the Julia job installs the candidate and runs
   `Pkg.test()` with the held-out suite.

Failure logs — including any expected values Test.jl prints — stay in this
private repository. Repository privacy replaces log redaction; no blind
evaluation harness is involved.

## Extending the suite

Each namespace has a directory under `test/full_api/`; add one file per public
object or tightly coupled helper family, then include it from `test/runtests.jl`.
The held-out goldens are a floor, not a ceiling: add randomized property checks
and real multi-call integration tests. `ci/check_full_api_plan.py` prevents the
source denominator from silently shrinking.
