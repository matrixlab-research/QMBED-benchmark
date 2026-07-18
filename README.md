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
                     key. Pre-filled from the frozen bundle's held-out goldens,
                     then hand-extended with properties + integration tests.
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

`test/runtests.jl` is a normal Julia test file — write whatever `@testset`,
`@test`, `@test_throws`, loops, and helpers you like. The held-out goldens are a
floor, not a ceiling: add randomized property checks and real multi-call
integration tests. The richer this file, the stronger the gate. Regenerate the
scaffold after a spec re-freeze with `minos build`, but hand-written tests are
yours to maintain.
