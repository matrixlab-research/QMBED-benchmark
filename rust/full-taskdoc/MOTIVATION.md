# quspin — full-package motivation

Status: design input for the complete Rust rewrite. This document is not the
frozen first-stage Minos task under `../taskdoc/` and does not change that
accepted contract.

Read this file together with `CONTRACT.md`, which defines the proposed public
boundary, and `TESTS.md`, which defines what must pass before a capability is
called complete.

## Why a full Rust package is needed

The current Rust crate proves that a small, Rust-native exact-diagonalization
core can execute the twelve paper-shaped workflows used by this repository.
That is a useful implementation milestone, but it is not equivalent to the
documented Python QuSpin 1.0.1 surface.

The full package should let a scientist:

1. describe a finite many-body Hilbert space and compatible symmetry sectors;
2. express static, time-dependent, and parameterized operators once;
3. choose stored or matrix-free execution without changing the physics;
4. solve spectra, evolution, Floquet, thermal-Lanczos, and measurement tasks;
5. move states and operators between sectors and composite bases;
6. reproduce the semantic results of the pinned Python reference on small
   systems and remain usable at paper-scale dimensions.

The goal is capability parity, not a transliteration of Python classes or
Julia multiple dispatch. Rust naming, ownership, iterators, traits, and error
handling should remain idiomatic where they do not change scientific meaning.

## Definition of complete

Completion has four independent layers:

1. **Surface coverage.** Every object in the frozen 64-object denominator has
   an explicit Rust mapping or an explicit, tested replacement.
2. **Semantic coverage.** Valid inputs give the same physical result, and
   invalid inputs fail in a stable documented category.
3. **Workflow coverage.** The public API composes into the twelve existing
   paper workflows without a private implementation hidden in the adapter.
4. **Scale usability.** Representative sparse, symmetry-reduced, and
   matrix-free cases avoid dense parent-space or full-Hilbert intermediates.

Passing the current 23-symbol task establishes a workflow core. It does not,
by itself, establish any of the four full-package layers above.

## Scientific capability inventory

### Bases and state representations

The package needs built-in spin, boson, spinless-fermion, spinful-fermion,
photon, tensor-product, and callback-defined bases. One-dimensional builders
cover common lattice sectors; general builders accept finite symmetry maps
without assuming a specific lattice or model.

State representation is a backend choice. A built-in fixed-width backend may
be fastest for small systems, but public semantics must not stop at `u128`.
Wide states are required for the Python `uint256`, `uint1024`, `uint4096`, and
`uint16384` denominator and for local actions on high site indices.

### Symmetry reduction

A symmetry sector is defined by finite maps, their periods, phases, and target
characters. Compatible commuting maps may be combined. The resulting basis
must expose deterministic representatives, normalization, projection, and
lifting, but physical results must not depend on the representative ordering.

The scale contract is as important as the mathematical contract: finite
permutation-like maps must be reduced through orbit metadata or another sparse
method. A valid large sector must not silently require a dense projected
symmetry matrix.

### Operators

Local operator strings and coupling lists form the reusable input language.
The same parsed terms must support:

- square Hamiltonians and rectangular sector-changing probes;
- dense, CSC, CSR, DIA, and matrix-free execution;
- real or complex coefficients;
- static, callable time-dependent, and parameterized components;
- bosonic factors, spin conventions, fermionic signs, and user callbacks.

Stored assembly and matrix-free action must share local transition semantics.
No workflow-specific assembler should define a second physical meaning for an
operator string.

### Solvers and dynamics

Dense diagonalization, selected Hermitian eigenpairs, shift-invert,
Lanczos decompositions, exponential action, general evolution, and Floquet
workflows all consume the same linear-operator boundary. Time-dependent
algorithms consume a corresponding time-dependent boundary instead of
rebuilding an unrelated operator abstraction.

Solver output includes residuals, convergence metadata, and the conventions
needed to interpret eigenvectors or invariant subspaces. Silent partial
success is not acceptable.

### Measurements and ensembles

The complete package includes expectation values, matrix elements,
fluctuations, partial traces, entanglement entropy, time-series observables,
diagonal ensembles, level statistics, subspace fidelity, state tracking, and
spectral response. Measurements must accept stored and matrix-free operators
where the underlying quantity only needs operator action.

### Persistence and interoperability

Operators need a versioned archive representation that preserves dense and
sparse values, shapes, and metadata. The first compatibility target is the
Python QuSpin ZIP/NPZ archive behavior already exercised by the Julia
verification suite. Loading untrusted archives must not execute code.

## First-principles design

### One narrow waist for state transitions

`Basis` is the minimal relation between a physical state and an index, plus
the result of applying one local operator. Built-in and user-defined bases use
this relation. Enumeration, sparse assembly, matrix-free action, and
cross-sector maps must not disagree about local physics.

### One narrow waist for linear maps

`LinearOperator` exposes shape and application. Dense, sparse, diagonal,
matrix-free, projected, exponential, Floquet, and Lindblad maps can implement
it. Algorithms should request additional capability only when mathematically
necessary, for example a reusable shifted solve for shift-invert.

### Time dependence is explicit

A time-dependent Hamiltonian is not a static operator with hidden mutable
global time. Its public boundary evaluates or applies the operator at a given
time. This makes parallel evaluation, caching, and reproducibility explicit.

### Costs are observable

Users must be able to tell whether an operation will enumerate a parent basis,
materialize a matrix, allocate a projector, or stay matrix-free. Public APIs
should not promise a storage format while constructing a dense intermediate.

### Capabilities, not model names, select algorithms

Optimized paths may depend on facts such as fixed particle number,
permutation-like symmetry, deterministic single-transition local action, or a
known matrix bandwidth. They must not depend on the workload being called XXZ,
SSH, PXP, or Hubbard.

## Compatibility with the current crate

The current public types remain the bootstrap layer:

- `Basis` and the five implemented basis builders;
- `OperatorTerm`, `OperatorBuilder`, `Operator`, and `LinearOperator`;
- `eigsh`, `evolve`, `Floquet`, `spectral_function`,
  `subspace_fidelity`, and `LindbladGenerator`.

The full design should extend these interfaces or provide mechanical migration
paths. It should not break the twelve verified workflows merely to make names
look closer to Python.

## Performance model

The following are design invariants, not optional optimizations:

- sparse assembly scales with generated transitions and stored nonzeros, not
  the square of the full Hilbert-space dimension;
- symmetry construction scales with orbit/projector data for finite maps;
- cross-sector action streams source transitions into target coordinates;
- matrix-free algorithms do not materialize the operator unless explicitly
  requested;
- repeated solves may reuse compatible factorizations and workspaces;
- batched state operations reuse typed buffers where practical;
- all scale claims are checked by allocation or peak-memory tests in addition
  to wall time.

## Non-goals

- Exact preservation of Python argument spelling, mutable attributes, or
  runtime type predicates.
- A single backend that simultaneously solves CPU, GPU, distributed, and
  out-of-core execution in the first complete release.
- Model-specific constructors that bypass the common basis/operator contract.
- Claiming complete-paper reproduction from paper-shaped workflow tests.
- Adding forward capabilities such as FQHE-specific bases or non-Abelian SU(2)
  sectors to the Python-1.0.1 parity denominator. They remain separate
  extensions.

## Delivery order

1. Remove current basis and solver restrictions while preserving the existing
   workflow API.
2. Complete general, composite, photon, projection, and wide-state bases.
3. Complete static/dynamic/parameterized operator semantics and persistence.
4. Complete solver, Floquet, block, Lanczos, measurement, and misc tools.
5. Activate Python differential tests for all mapped objects, then add scale
   gates that are independent of Python performance.

Each stage must include public behavior tests and private verification before
the corresponding row in the full contract is marked complete.
