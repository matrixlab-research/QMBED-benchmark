# Fifty research applications of exact diagonalization

This document is a workflow-led research and benchmark horizon for QMBED. It
does **not** claim that the repository already executes fifty new benchmarks.
The current executable suite remains the 21 deterministic integration scenarios
and twelve paper-shaped timed workflows documented elsewhere.

The machine-readable source of truth is
[`benchmark/catalog/lkm_ed_applications.json`](../benchmark/catalog/lkm_ed_applications.json).
Each catalog entry states a candidate medium-scale problem, primary observable,
numerical validation, current QMBED boundary, automatic-differentiation (AD)
role, and the capabilities required before promotion into an executable
benchmark.

## Discovery method and interpretation

The snapshot was collected on 2026-07-27 from the
[LKM Workflow Families index](https://lkm.bohrium.com/web/zh/workflows).
The search used 94 query variants spanning exact diagonalization, Lanczos,
many-body models, dynamics, response functions, topology, open systems,
quantum control, inverse problems, and automatic differentiation. The queries
returned 603 distinct workflow families. Forty-six high-signal families were
then checked against their complete definition, step skeleton, inputs/outputs,
technique choices, candidate tasks, numerical validation protocol, and paper
count.

Four additional LKM claim searches and four reasoning-chain searches were used
for differentiable ED. They exposed application chains that are not named as
standalone workflow families, including dominant-eigensolver reverse mode,
Hamiltonian learning, inverse entanglement design, differentiable sparse
dynamics, and quantum optimal control.

LKM cluster IDs below record **discovery provenance**, not scientific
citations. Representative DOI sources are retained in the JSON catalog. A
retrieval score was used only for ranking LKM results; it was not interpreted
as scientific confidence.

The coverage labels are intentionally strict:

- **Direct**: QMBED can express the model, solver path, and main observable.
- **Partial**: core ED kernels exist, but a domain builder, sweep, or
  post-processing layer is missing.
- **Beyond**: at least one load-bearing scientific or differentiation
  capability is absent.

The resulting denominator is 14 direct, 15 partial, and 21 beyond-QMBED
applications. Thirteen applications require AD, rather than merely benefiting
from it.

## 1. Equilibrium states and phase diagrams

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-001 | Symmetry-resolved Heisenberg/XXZ spectra | 3831, 4905 | Direct | None | `L=20`, `Nup=10`, `k=0`, four eigenpairs |
| ED-002 | Frustrated `J1-J2` spin gap and phase boundary | 23806, 3831 | Direct | Helpful | `L=18`, two spin sectors, nine couplings |
| ED-003 | TFIM criticality from gap and fidelity | 30690, 23806 | Direct | Helpful | `L=18`, seven fields, tracked two-state subspace |
| ED-004 | Magnetization plateaus | 23177, 23806 | Partial | Helpful | `L=16`, all `Sz` sectors, 21 fields |
| ED-005 | Multi-orbital Hubbard phase diagram | 8704, 3127 | Partial | Helpful | Four-site/two-orbital cluster, `5x5` interaction grid |

These cases move from a single diagonalization to cross-sector envelopes,
tracked eigen-subspaces, and parameterized phase maps. The missing layer is
mostly orchestration; the exception is the ergonomic construction of general
multi-orbital local spaces.

## 2. Dynamics and transport

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-006 | Bias-release Hubbard transport | 32600 | Direct | Helpful | Spinful `L=10`, dimension 63504, 101 times |
| ED-007 | Floquet heating and quasienergy analysis | 1069 | Direct | Helpful | Driven Ising `L=10`, full unitary |
| ED-008 | Strong-coupling doublon dynamics | 31254 | Direct | Helpful | Two particles on `L=14`, six `U/t` values |
| ED-009 | Many-body scar states and revivals | 22708 | Direct | Helpful | Periodic PXP `L=24`, dimension 103682 |
| ED-010 | Perfect state transfer in a spin chain | 30037 | Direct | Helpful | `L=20` one-excitation sector with disorder |

The forward simulations are already close to the QMBED narrow waist:
parameterized operators, sparse action, Krylov propagation, and observables.
The frontier is the backward pass through propagation and reusable ensemble
execution.

## 3. Response functions and spectroscopy

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-011 | Particle-addition/removal spectral function | 1848 | Direct | Helpful | Eight-site Hubbard, cross-sector Lanczos |
| ED-012 | Dynamical spin structure factor | 7817, 3831 | Direct | Helpful | `L=18`, momentum source, 150 Lanczos steps |
| ED-013 | RIXS on Hubbard/Cu-O clusters | 5881, 3127 | Partial | Helpful | Ground/core/final sectors, 2D energy grid |
| ED-014 | Optical conductivity and Drude weight | 968 | Partial | Helpful | Half-filled `L=18`, current response |
| ED-015 | Impurity Green function and Kondo crossover | 30508, 1848 | Partial | Helpful | One impurity plus six spinful bath orbitals |

The universal kernel is a source state plus a resolvent or Krylov spectral
measure. RIXS, conductivity, and impurity workflows should therefore add
composable domain operators and sum-rule checks, not workflow-specific
eigensolvers.

## 4. Topology and entanglement

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-016 | Interacting SSH transition | 19000, 16282 | Direct | Helpful | `L=18`, seven dimerizations, three interactions |
| ED-017 | FQHE spectra and trial-state overlaps | 12454 | Beyond | Helpful | Eight electrons on a sphere |
| ED-018 | Moiré fractional Chern insulator | 26341 | Beyond | Helpful | Projected `6x4` Chern band with twists |
| ED-019 | Majorana parity splitting | 5931, 21745 | Partial | Helpful | `L=16` even/odd parity sectors |
| ED-020 | Adiabatic braiding and Berry matrices | 27880, 5931 | Beyond | Helpful | Six to eight Majoranas, 100-point path |

FQHE/FCI and braiding are not ordinary spinless-fermion benchmarks. They require
magnetic or projected-band matrix elements, two-dimensional momentum sectors,
degenerate-manifold tracking, and gauge-covariant Berry/Wilson-loop
observables. The catalog keeps those boundaries explicit.

## 5. Open systems and quantum optics

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-021 | Collective emission and subradiance | 14428 | Direct | Helpful | Eight emitters, matrix-free Liouvillian |
| ED-022 | Lindblad steady state and gap | 231, 195 | Direct | Helpful | Dissipative `L=8`, Liouville dimension 65536 |
| ED-023 | Quantum-trajectory unraveling | 8037, 231 | Beyond | Helpful | `L=10`, 1000 seeded trajectories |
| ED-024 | Cavity-QED photon blockade | 10264, 195 | Partial | Helpful | Two emitters, photon cutoff 16 |
| ED-025 | Dicke superradiant transition | 30105 | Partial | Helpful | Four spins, photon cutoffs 16/24/32 |

QMBED has a matrix-free Lindblad foundation, but trajectory sampling,
trace-constrained steady-state solves, correlation helpers, and systematic
local-cutoff convergence remain separate capabilities.

## 6. Cross-domain ED

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-026 | Few-electron quantum-dot CI | 2880 | Partial | Helpful | Twelve Fock-Darwin orbitals, four electrons |
| ED-027 | Nuclear shell-model spectroscopy | 4003, 4905 | Beyond | Helpful | `sd` shell, dimension about `1e5` |
| ED-028 | Molecular full-CI | 6178, 463 | Beyond | Helpful | `CAS(10e,10o)`, six states |
| ED-029 | Truncated-Fock-space quantum field theory | 21731 | Beyond | Helpful | Dimension about `1e5`, four cutoffs |
| ED-030 | Self-consistent BdG vortex/LDOS | 22187, 19234, 24103 | Beyond | Helpful | `24x24` lattice, up to 50 iterations |

These examples share sparse eigensolvers but not the same scientific basis
semantics. They are strong candidates for domain packages above a small QMBED
operator/eigensolver core, rather than reasons to enlarge the core with five
unrelated special paths.

## 7. Thermodynamics, localization, and certification

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-031 | Finite-temperature susceptibility fitting | 11058, 3831 | Partial | Helpful | `L=16` full spectrum or `L=24` FTLM |
| ED-032 | Thermal typicality and thermalization | 12030 | Partial | Helpful | `L=18` with exact `L=12` oracle |
| ED-033 | Anderson localization and multifractality | 19021, 29318 | Direct | Helpful | `20^3` lattice, 20 interior eigenpairs |
| ED-034 | MBL level statistics | 19021, 3831 | Partial | Helpful | `L=18`, 50 disorder samples |
| ED-035 | ED certification of variational ansatzes | 12290, 3831 | Partial | Helpful | `L=20`, batch of 100 trial states |

The dominant missing capability is not a new Hamiltonian assembler. It is
reproducible stochastic/ensemble orchestration: FTLM or typicality vectors,
parallel disorder samples, confidence intervals, and resource-aware
shift-invert.

## 8. Metrology and control

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-036 | Entropy and mixed-state negativity | 1989 | Partial | Helpful | Pure `L=18`; mixed `L=8` |
| ED-037 | Many-body quantum Fisher information | 18404, 30690 | Partial | Helpful | `L=16`, pure and thermal QFI |
| ED-038 | Closed-system pulse optimization | 584 | Beyond | **Essential** | Six transmons or 16 sparse spins |
| ED-039 | Dissipative state preparation | 195, 231 | Beyond | **Essential** | Two qubits plus cavity, 100 controls |
| ED-040 | Many-body sensing optimization | 18404, 195 | Beyond | **Essential** | `L=16`, 100 controls |

The benchmark contract for control must time both value and gradient, check
directional derivatives, and report memory/checkpoint costs. A forward-only
simulation benchmark is not evidence that an optimization workflow is usable.

## 9. Differentiable inverse problems

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-041 | AD fidelity susceptibility | 30690, 3831 | Beyond | **Essential** | TFIM `L=20`, dominant eigenpair only |
| ED-042 | Learn spin Hamiltonian from thermodynamics | 11058, 3831 | Beyond | **Essential** | `L=12`, five couplings, 50 temperatures |
| ED-043 | Inverse design of entangled ground states | 1989, 3831 | Beyond | **Essential** | `L=16`, 20-40 couplings |
| ED-044 | Transmon-resonator inverse spectroscopy | 584, 11058 | Beyond | **Essential** | 50-200 dimensional circuit Hamiltonian |
| ED-045 | Hamiltonian reconstruction from time traces | 32600 | Beyond | **Essential** | Eight sites, 20 coefficients, 100 times |

These applications convert ED from a terminal solver into a differentiable
model inside an optimizer. Their validation therefore includes synthetic
parameter recovery and held-out prediction, not only a gradient check.

## 10. Differentiable numerical kernels

| ID | Application | LKM family | Coverage | AD | Candidate medium-scale shape |
|---|---|---:|---|---|---|
| ED-046 | Reverse-mode dominant sparse eigensolver | 3831 | Beyond | **Essential** | Dimension about `1e6`, `k=1-4` |
| ED-047 | Degenerate-projector differentiation | 26341, 30690 | Beyond | **Essential** | Two-/fourfold manifolds in dimension `1e5` |
| ED-048 | Differentiable sparse Krylov evolution | 32600, 584 | Beyond | **Essential** | `L=18`, 200 times, 100-200 parameters |
| ED-049 | Differentiable impurity/NRG iteration | 30508 | Beyond | **Essential** | Impurity plus eight bath sites |
| ED-050 | Differentiable CI/few-body basis optimization | 6178, 463 | Beyond | **Essential** | `CAS(8e,8o)` or 50-100 basis parameters |

These five kernels define the reusable AD waist:

1. a parameterized `LinearOperator` with matrix-vector, adjoint action, and
   parameter VJP/JVP;
2. dominant-eigenpair reverse mode through projected linear solves, without a
   full spectrum;
3. degenerate-subspace outputs represented by projectors rather than arbitrary
   eigenvector gauges;
4. checkpointed reverse mode for Krylov/exponential evolution; and
5. stable generalized-eigenproblem and truncation semantics.

The first implementation milestone should cover ED-046 through ED-048 because
they unlock multiple applications in sections 8 and 9 without introducing a
workflow-specific path. Domain-specific CI, NRG, FQHE/FCI, nuclear, and BdG
support can then compose over the same boundary.

## Representative AD evidence returned by LKM

- [Automatic differentiation of dominant eigensolver and its applications in
  quantum physics](https://doi.org/10.1103/PhysRevB.101.245139) motivates
  low-rank projected backward solves and higher-order derivatives.
- [Inverse Hamiltonian design of highly entangled quantum
  systems](https://doi.org/10.1103/PhysRevResearch.6.033080) exposes
  degeneracy and eigenvector-gauge failure modes in differentiable ED.
- [Learning Effective Spin Hamiltonian of Quantum
  Magnet](https://doi.org/10.48550/arXiv.2011.12282) uses a differentiable ED
  thermodynamics solver for parameter inference.
- [Achieving fast high-fidelity optimal control of many-body quantum
  dynamics](https://doi.org/10.48550/arXiv.2008.06076) requires gradients that
  are exact for the discrete propagation used in the forward pass.
- [Automatic differentiable numerical
  renormalization group](https://doi.org/10.1103/PhysRevResearch.4.013227)
  shows why iterative truncation needs explicit differentiable semantics.
- [Differentiable quantum chemistry with PySCF](https://doi.org/10.48550/arXiv.2207.13836)
  supplies generalized-eigenproblem and degeneracy cases for the CI frontier.

## Promotion rule

An entry becomes an executable benchmark only when it has:

1. deterministic inputs and a pinned scientific/method source;
2. a direct numerical oracle, sum rule, analytic limit, or independent
   implementation;
3. a medium-size resource envelope suitable for the selected CI/HPC tier;
4. explicit forward accuracy and, when applicable, gradient accuracy;
5. a language-neutral observation schema; and
6. no workflow-name dispatch or model-specific solver shortcut.

This keeps the catalog useful as an expansion map without weakening the
meaning of the existing green benchmark suite.
