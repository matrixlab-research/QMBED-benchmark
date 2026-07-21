# LKM-derived workflow validation

This suite separates fast package tests from held-out scientific validation.
The split is based on numerical role, not on model names:

- **QuSpin.jl CI tests** use tiny systems with analytic or independently
  constructed answers. They test public API contracts on every push and pull
  request.
- **Private validation** uses larger Hilbert spaces, parameter scans,
  truncation checks, cross-sector observables, and paper-shaped invariants.
  These tests live only in this repository and are timed by the verification
  workflow.

LKM supplied workflow and paper discovery. It is not a runtime dependency and
its internal family identifiers are not citations. The publications below are
the scientific/method sources used to select the validation shapes; the tests
are deterministic reduced instances, not paper reproductions.

## Public CI test set

| Capability | Public test | Fast invariant |
|---|---|---|
| Locally constrained state generation | `test/unit/workflow_analysis.jl` | six-site periodic blockade space has 18 states |
| Degenerate-subspace fidelity/tracking | same | invariant under an arbitrary unitary rotation inside the subspace |
| Dynamical correlator | same | two-level result is `exp(-im*gap*t)` |
| Lehmann/resolvent/Krylov spectrum | same | all methods reproduce one Lorentzian transition |
| Matrix-free Lindblad generator | same | analytic amplitude-damping population, trace, Hermiticity, positivity |

The small performance suite also records paired `L=20` PXP state-generation
rows for prefix pruning and full-space predicate filtering. These timings are
observational and do not introduce a noisy CI threshold.

## Private validation set

| Validation workflow | Size/path | Verified scientific or numerical invariant | Source |
|---|---|---|---|
| TFIM fidelity scan | `L=10`, six fields | adjacent fidelity minimum and finite-size gap response near the transition | LKM family 30690, fidelity-based transition detection |
| Kitaev parity-sector gap | `L=8`, even/odd Fock sectors | exponentially small topological parity splitting versus a finite trivial splitting | [Sau and Das Sarma](https://doi.org/10.48550/arXiv.1111.6600) |
| PXP scar revival | periodic `L=14`, constrained dimension 843 | norm conservation and revival of the staggered state | [Pal et al.](https://doi.org/10.48550/arXiv.2411.02500) |
| Bose-Hubbard Mott quench | `L=7`, `Nb=7`, `sps=3` | number/norm conservation and decay of the Mott return probability | [Queisser et al.](https://doi.org/10.1103/PhysRevA.89.033616) |
| Loschmidt return | TFIM `L=10` quench | nontrivial return-rate growth after crossing the critical region | [Mueller et al.](https://doi.org/10.48550/arXiv.2210.03089) |
| Dicke cutoff convergence | three spins, photon cutoffs 4 and 6 | ground energy and photon occupation converge with truncation | [Lewis-Swan et al.](https://doi.org/10.1103/PhysRevResearch.3.L022020) |
| Spinful-Hubbard transport | half-filled `L=6`, dimension 400 | bias-induced current and relaxation of density imbalance | LKM family 32600, time-dependent lattice transport |
| Doublon dynamics | two particles on `L=6` | large-`U` doublon survival exceeds the free result | LKM family 31254, strong-coupling doublon dynamics |
| Dynamical spin structure factor | Heisenberg `L=10`, dimension 252 | positivity, spectral sum rule, and finite-frequency peak | [Shimokawa et al.](https://doi.org/10.48550/arXiv.2206.13064) |
| Particle-addition spectrum | spinless `L=8`, `N=3 -> 4` | cross-sector spectral weight equals the operator sum rule | [Pichler et al.](https://doi.org/10.48550/arXiv.2410.07319) |
| OTOC growth | TFIM `L=7` | initially commuting distant operators develop a nonzero commutator | [Li et al.](https://doi.org/10.1103/PhysRevResearch.2.043399) |
| Collective Lindblad decay | four emitters, density matrix dimension 16 | bright-state population follows `exp(-N*gamma*t)` while trace and Hermiticity are preserved | collective-emission workflow; matrix-free Lindblad extension |

The RIXS and Anderson-impurity workflows are represented at the shared
response-kernel level by the same-sector and cross-sector spectral tests.
Full Kramers-Heisenberg intermediate-state modeling and distributional bath
fitting remain application-layer workflows, not claims of direct package
coverage.
