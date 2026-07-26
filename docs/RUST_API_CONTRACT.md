# Rust QMBED API contract

Status: implemented paper-workflow core and bootstrap design for the full Rust
package.

The accepted first-stage Minos contract is frozen under `rust/taskdoc/`. The
complete Python-1.0.1 capability design is now specified separately in:

- [`rust/full-taskdoc/MOTIVATION.md`](../rust/full-taskdoc/MOTIVATION.md);
- [`rust/full-taskdoc/CONTRACT.md`](../rust/full-taskdoc/CONTRACT.md);
- [`rust/full-taskdoc/TESTS.md`](../rust/full-taskdoc/TESTS.md).

Those documents are the implementation-preparation boundary for the full
package. This file remains the shorter architectural overview and migration
map for the original paper-workflow core.

## Design boundary

The Rust package should not imitate Python or Julia syntax. It should preserve
the scientific semantics that the verification corpus actually exercises:

1. define a Hilbert-space basis;
2. parse reusable local operator terms;
3. either assemble a requested matrix format or expose matrix-free action;
4. run eigensolvers, time evolution, spectra, and observables against the same
   linear-operator interface;
5. report explicit errors for invalid symmetries, particle sectors, dimensions,
   and unsupported storage.

The narrow waist is `Basis` plus `LinearOperator`. CSC is one backend, not the
definition of a Hamiltonian. This keeps the PXP constrained-basis work useful
for spin, boson, fermion, symmetry-resolved, and future GPU/MPI implementations.

## Proposed public modules

```text
qmbed/
  basis        Basis, SpinBasis1D, BosonBasis1D, SpinlessFermionBasis1D,
               SpinfulFermionBasis1D, PhotonBasis, TensorBasis, UserBasis
  operator     OperatorSpec, OperatorBuilder, Hamiltonian, MatrixFormat,
               LinearOperator, operator_matrix, sector-changing operators
  solve        eigsh, lanczos, expm_multiply, evolve
  dynamics     Floquet, dynamical_correlator, spectral_function
  measure      expectation, entanglement_entropy, diagonal_ensemble
  workflow     subspace_fidelity, track_eigenspaces, LindbladGenerator
  error        QmbedError, Result
```

## Core interfaces

The exact ownership details may change during implementation, but these
semantics should remain stable.

```rust
pub trait Basis: Send + Sync {
    type State: Copy + Eq + Send + Sync;

    fn len(&self) -> usize;
    fn state(&self, index: usize) -> Result<Self::State>;
    fn index(&self, state: Self::State) -> Result<usize>;

    // Applies one parsed local operator to one basis state. The assembler and
    // matrix-free path share this primitive.
    fn apply_local(
        &self,
        state: Self::State,
        operator: &ParsedOperator,
        sites: &[usize],
    ) -> Result<Option<(Self::State, Complex64)>>;
}

pub trait LinearOperator: Send + Sync {
    fn shape(&self) -> (usize, usize);
    fn apply(&self, input: &[Complex64], output: &mut [Complex64]) -> Result<()>;
}
```

`shape()` is `(rows, columns)`, so this same interface also represents a
sector-changing probe. Hamiltonians and time generators are square; a
particle-addition operator may map an N-particle source basis to an
(N+1)-particle target basis. Algorithms validate the shape they require.

Algorithms accept `&impl LinearOperator`; therefore the same `eigsh`, Krylov
evolution, spectrum, and Lindblad code works with dense, CSC, CSR, DIA, or a
matrix-free operator. Materialization is explicit:

```rust
let h = OperatorBuilder::on(&basis)
    .terms(terms)
    .checks(AssemblyChecks::all())
    .build(MatrixFormat::Csc)?;

let eigenpairs = eigsh(
    &h,
    EigshOptions::smallest_algebraic(6)
        .with_tolerance(1e-9)
        .with_krylov_dimension(32),
)?;
```

## Basis construction

Use typed builders for built-in bases and a callback-based `UserBasis` for
constraints. Site indices are zero-based in Rust. State enumeration order is
not part of the public contract; tests use physical states or
`basis.index(state)` rather than hard-coded vector positions.

```rust
let spin = SpinBasis1D::builder(18)
    .magnetization(9)
    .momentum(0)
    .pauli(false)
    .build()?;

let pxp = UserBasis::builder(24)
    .state_filter(periodic_blockade)
    .operator('x', spin_flip)
    .build()?;
```

`UserBasis` must feed the same universal assembly path as built-in bases. A
PXP-only CSC constructor would pass one benchmark but fail the architecture.

## Operator terms and assembly

```rust
let terms = vec![
    OperatorSpec::new("+-", hopping_forward)?,
    OperatorSpec::new("-+", hopping_backward)?,
    OperatorSpec::new("zz", interactions)?,
];
```

- Operator strings are parsed and validated once.
- Couplings carry a complex coefficient and a fixed number of sites.
- Fermionic sign handling belongs to the basis/operator application layer.
- Assembly streams triplets into a format-specific sink; it must not create a
  dense matrix before CSC/CSR materialization.
- Hermiticity, particle conservation, and symmetry checks are individually
  configurable and return structured errors.

## Solvers and workflows

```rust
pub fn eigsh(op: &impl LinearOperator, options: EigshOptions)
    -> Result<Eigensystem>;

pub fn evolve(
    op: &impl LinearOperator,
    initial: &[Complex64],
    options: EvolutionOptions,
) -> Result<StateTrajectory>;

pub fn spectral_function(
    h: &impl LinearOperator,
    source: &[Complex64],
    probe: &impl LinearOperator,
    options: SpectrumOptions,
) -> Result<Vec<f64>>;

pub fn subspace_fidelity(left: &Subspace, right: &Subspace) -> Result<f64>;
```

Important semantic requirements:

- `eigsh` exposes target selection, optional shift-invert, deterministic seed,
  tolerance, maximum iterations, and residuals.
- `evolve` supports a single time and a time grid without first densifying.
- subspace fidelity is invariant under unitary rotations inside a degenerate
  subspace.
- spectral functions support same-sector and cross-sector sources; the probe
  shape is `(target_dimension, source_dimension)`.
- `LindbladGenerator` implements `LinearOperator` over vectorized density
  matrices, so open-system evolution can remain matrix-free.

## Error model

All fallible public operations return `Result<T, QmbedError>`. The minimum
stable categories are:

- invalid operator string or coupling arity;
- state/sector not present in the basis;
- incompatible symmetry or particle-number sector;
- non-Hermitian input when a Hermitian solver is requested;
- dimension mismatch;
- unsupported storage/solver combination;
- solver non-convergence with iteration count and final residual.

Panics are reserved for internal invariant violations, not user input.

## Verification adapter

The public verification crate under `rust/` owns two traits:

- `QmbedApi`: a low-level adapter over basis construction, Hamiltonian
  assembly, solvers, dynamics, spectra, and Lindblad evolution;
- `WorkflowBackend`: an end-to-end adapter that returns named physical metrics
  for the twelve paper workflows.

The adapter remains external to the candidate crate. This prevents the public
package from depending on its verification repository and avoids forcing
public types to implement a foreign verification trait.

The original Python-to-Julia denominator remains the source of truth. The
Rust-side [`full_api_contract.json`](../rust/full_api_contract.json) maps all
64 objects, 282 non-constructor methods, and 180 attributes into `basis`,
`operator`, `solve`, `dynamics`, `measure`, and `workflow`. CI rejects changes
that shrink the source denominator, leave a namespace unmapped, or reference a
missing Rust contract test. Numerical activation advances separately through
`adapter_compiles`, `oracle_passes`, `workflow_passes`, and
`paper_benchmark_recorded`.

## Migrated benchmark denominator

The Rust catalog preserves the exact case IDs and parameters already used for
Python and Julia:

| Group | Workflows |
|---|---|
| Original six | MBL shift-invert; XXZ Lanczos quench; Floquet full unitary; spinful Hubbard; interacting SSH; translation-sector XXZ |
| LKM extension | TFIM subspace fidelity; PXP revival; Bose-Hubbard Mott quench; Hubbard current quench; CoNb2O6 response; particle-addition spectrum |

Each case declares required capabilities and named physical invariants. Timing
is emitted only after the invariants pass. The CSV schema matches the current
renderer and adds `language=rust`, so Python/Julia/Rust bars can share one
report without special cases.

## Julia-to-Rust migration map

| Julia verification call | Proposed Rust surface |
|---|---|
| `SpinBasis1D(...)` | `SpinBasis1D::builder(...).build()?` |
| `UserBasis(...; states=...)` | `UserBasis::builder(...).state_filter(...).operator(...).build()?` |
| `OperatorSpec(opstr, couplings)` | `OperatorSpec::new(opstr, couplings)?` |
| `Hamiltonian(basis, terms; static_fmt=:csc)` | `OperatorBuilder::on(&basis).terms(terms).build(MatrixFormat::Csc)?` |
| `H * vector` | `LinearOperator::apply(&h, input, output)?` |
| `eigsh(H; ...)` | `eigsh(&h, EigshOptions::new(k, target).with_tolerance(tol))?` |
| `lanczos_full` + `expm_lanczos` | `evolve(&h, &psi0, EvolutionOptions::new(times))?` |
| `track_eigenspaces` | `track_eigenspaces(&subspaces)?` |
| `spectral_function` | `spectral_function(&h, &source, &probe, options)?` |
| `LindbladGenerator` | matrix-free `LindbladGenerator` implementing `LinearOperator` |

## Staged implementation

1. Implement spin bases, `OperatorSpec`, universal COO-to-CSC assembly,
   `LinearOperator`, matvec, and low-energy `eigsh`.
2. Add fermion/boson/user bases without changing solver interfaces.
3. Add Krylov evolution, Floquet, subspace tracking, and spectra.
4. Add matrix-free Lindblad and advanced symmetry/storage backends.

Every stage should activate the corresponding independent tests and paper
workflows; API presence alone is not completion.
