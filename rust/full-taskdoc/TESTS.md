# quspin — full-package tests

Status: acceptance design for the complete Rust package. These tests extend,
rather than replace, the frozen first-stage examples in `../taskdoc/TESTS.md`.

## What the test program must establish

The suite must distinguish four questions:

1. Does the public Rust surface exist?
2. Does it implement the intended mathematics?
3. Does it compose into real scientific workflows?
4. Does it remain usable at representative sparse and symmetry-reduced scale?

A green answer to one question must not be reported as a green answer to the
others.

## Public/private split

Public crate tests contain deterministic examples, algebraic properties,
regressions, and small scale checks. The private verification repository owns
independent parameters, held-out seeds, Python oracle comparison, and
paper-workflow release gates.

The private adapter may translate public Rust outputs into observations. It
must not implement basis enumeration, Hamiltonian assembly, eigensolvers,
evolution, or measurements on behalf of the candidate.

## Gate 0 — denominator and documentation

- The frozen Python denominator remains 64 objects, 282 non-constructor
  methods, and 180 attributes unless an explicitly reviewed source update is
  made.
- Every Python object maps to a Rust symbol or a documented many-to-one
  replacement in `CONTRACT.md`.
- Every mapped object names at least one public test and one private semantic
  family.
- The first-stage 23-symbol contract remains separately reproducible.
- CI rejects removal of a mapping, test family, or paper-workflow requirement.

## Gate 1 — basis semantics

| Capability | Public properties | Independent verification |
|---|---|---|
| State/index bijection | round trip every state in small sectors | held-out sectors and wide-state samples |
| Spin bases | combinatorial dimensions, spin/Pauli conventions, higher-spin ladder algebra | Python dimensions, projectors and operator matrices |
| Boson bases | cutoff and particle-number dimensions, ladder factors | Python basis and Hamiltonian observations |
| Fermion bases | canonical signs, species ordering, particle-hole maps | independent Clifford/sign properties and Python matrices |
| General symmetries | map period, character, projector orthonormality, commuting-map intersection | held-out 1D/2D maps and Python projectors |
| Tensor basis | state/index factorization and multi-factor Kronecker action | independent explicit product-space oracle |
| Photon basis | oscillator ladder action and total-excitation sectors | Python dimensions, matrices and entropy |
| User basis | deterministic filters, callbacks, map sectors, deferred construction | held-out callbacks and direct transition tables |
| Wide states | high-site bit action, conversion, shift and mask identities | randomized round trips for every fixed width |
| Projection | lift/project identities and particle-conserving output modes | independent sparse projector multiplication |
| Partial trace/entropy | trace, positivity and product/maximally-entangled limits | independent dense reduced-density oracle |
| Sector shift | stored and matrix-free agreement for vector and batch inputs | explicit target-projector times full action times source-projector |

Invalid sites, sectors, map periods, incompatible characters, state overflows,
and unsupported local operators must return stable errors and never panic.

## Gate 2 — operator semantics

| Capability | Public properties | Independent verification |
|---|---|---|
| Parsed terms | arity and site validation; coefficient linearity | held-out operator strings and complex coefficients |
| Stored formats | Dense/CSC/CSR/DIA values agree; sparse indices canonical | Python matrices plus format-structure assertions |
| Matrix-free action | agrees with stored action for states and batches | independently generated transition sums |
| Dynamic Hamiltonian | evaluation at multiple times; transform/algebra consistency | original and transformed Python operators |
| Parameterized operator | missing/default parameters and component linearity | explicit weighted component oracle |
| Cross-sector maps | correct rectangular shape, signs and normalization | independent projector formula |
| Operator algebra | sum, product, powers, commutator and anticommutator identities | dense small-system oracle |
| Observables | expectation, matrix element and variance identities | direct dense vector contractions |
| Conversions | values preserved across supported formats and element types | round-trip structure and value checks |
| Archives | Rust round trip and bidirectional Python interoperability | held-out dense/sparse NPZ fixtures |

Sparse constructors are instrumented to reject a dense full-Hilbert
intermediate. Matrix-free constructors are instrumented to reject any complete
matrix materialization.

## Gate 3 — solvers and evolution

### Dense and partial eigensolvers

- Real-symmetric and complex-Hermitian cases satisfy residual and
  orthonormality tolerances.
- Algebraic, magnitude, both-end, and shift targets return the requested
  invariant subspace.
- Results are compared as invariant subspaces near degeneracy, not by arbitrary
  eigenvector phase or rotation.
- Matrix-backed and matrix-free paths agree.
- Initial vectors, deterministic seeds, iteration limits, and
  `return_eigenvectors=false` behavior are exercised.
- Forced non-convergence returns structured iteration/residual metadata.

### Lanczos and exponential action

- `lanczos_full` reconstructs the Krylov projection and reports an orthonormal
  basis within tolerance.
- `lanczos_iter` yields the same sequence without storing all vectors.
- `expm_lanczos` and general `expm_multiply` agree with dense exponential
  action on small systems.
- Reusable plans/workspaces produce the same result as fresh calls.
- FTLM and LTLM estimators agree with exact small-system thermal traces within
  their documented statistical tolerance.

### Evolution

- Time zero returns the initial state.
- Composition over adjacent intervals agrees with direct evolution.
- Hermitian Hamiltonian evolution preserves norm.
- Liouville-von Neumann evolution preserves trace and Hermiticity.
- Batched evolution equals independent column-wise evolution.
- Static and callable time-dependent generators are both checked.
- Collected and iterator output contain the same times and states.

## Gate 4 — Floquet, block, measurement, and workflow tools

| Capability | Required checks |
|---|---|
| Floquet time vector | cycle/time coordinates, endpoint conventions and monotonicity |
| Floquet full unitary | unitarity, quasienergy convention and eigenpair residuals |
| Floquet action | repeated period application agrees with explicit small-system unitary |
| Block tools | block spectra and dynamics reconstruct full-space results |
| Observable time series | each recorded value equals an independent call at the same state/time |
| Diagonal ensemble | normalization, energy-window selection and fluctuation formulas |
| Level statistics | filtering/degeneracy conventions and analytic toy spectra |
| Spectral response | Lehmann small-system oracle and same/cross-sector probes |
| Subspace fidelity | gauge invariance under unitary rotations inside each subspace |
| State tracking | permutation/phase invariance and explicit ambiguity diagnostics |
| Lindblad generator | trace and Hermiticity preservation; explicit-small-system agreement |

## Gate 5 — paper workflows

All twelve existing cases remain required:

1. MBL interior spectrum with shift-invert;
2. XXZ quench with Krylov evolution;
3. Floquet heating with full-unitary checks;
4. spinful Hubbard spectrum;
5. interacting SSH spectrum;
6. translation-sector XXZ spectrum;
7. TFIM degenerate-subspace fidelity scan;
8. constrained PXP revival;
9. Bose-Hubbard Mott quench;
10. spinful-Hubbard current quench;
11. CoNb2O6 dynamical response;
12. particle-addition cross-sector spectrum.

Each case declares required public capabilities, deterministic parameters,
named physical invariants, and its direct/assisted/proxy relationship to the
source literature. Timing begins only after the invariants pass.

The workflow gate is necessary but not sufficient for full-package parity:
the 64-object denominator includes capabilities that these twelve workflows do
not exercise.

## Gate 6 — scale and allocation

The following cases must be checked across multiple sizes, not one selected
benchmark point:

- fixed-particle spin enumeration;
- one-dimensional translation/parity construction;
- two-map general-lattice symmetry intersection;
- spin, boson, spinless, spinful, and user-basis CSC assembly;
- cross-sector particle and symmetry shifts;
- matrix-free action and Krylov evolution;
- sparse shift-invert with repeated compatible solves;
- batched evolution and entropy;
- wide-state local action above the `u128` range.

Acceptance is expressed structurally before it is expressed as speed:

- no full-square operator in sparse or cross-sector transition paths;
- no dense projected-symmetry fallback for permutation-like finite maps;
- no parent-space dense output when direct target accumulation is possible;
- no operator materialization in matrix-free solver paths;
- allocations and peak memory follow transition/projector nonzeros rather than
  full-space dimension squared.

Wall-time comparisons use same-runner Python/Rust measurements, warm-up plus
multiple samples, fixed thread counts, raw sample retention, and physical
preflights. Performance remains observational until runner variance supports a
stable regression threshold.

## Activation states

Every capability family moves independently through:

```text
planned
→ public_surface
→ public_properties
→ private_oracle
→ workflow_covered
→ scale_covered
→ complete
```

Namespace-level workflow success must not advance every object in that
namespace. Status is recorded per mapped object or tightly coupled semantic
family.

## Initial implementation batches

### Batch A — remove restrictions in the existing core

- higher spin, parity, and fermion momentum sectors;
- complex-Hermitian small eigensystems and the remaining spectrum targets;
- time-dependent and batched evolution;
- transition-driven cross-sector action.

### Batch B — complete bases

- four general bases and finite symmetry maps;
- tensor and photon bases;
- wide state backends;
- projection, partial trace, entropy, and advanced UserBasis lifecycle.

### Batch C — complete operators

- dynamic Hamiltonian and parameterized `QuantumOperator`;
- `ExpOp`, algebra, observables, conversions, and archives;
- reusable exponential and shifted-solve plans.

### Batch D — complete tools

- full Floquet, block tools, public Lanczos, FTLM/LTLM;
- measurements, diagonal ensembles, misc statistics and conversion helpers;
- state tracking and workflow extensions.

No batch is complete until its public tests, private verification, and relevant
scale checks are all active.
