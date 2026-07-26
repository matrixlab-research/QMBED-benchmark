# quspin — full-package contract

Status: proposed public boundary for the complete Rust package. The accepted
first-stage contract in `../taskdoc/CONTRACT.md` remains frozen and authoritative
for the existing 23-symbol gate.

## Authority and scope

The source denominator is the pinned Python QuSpin 1.0.1 public surface in
`../../test/full_api_migration_plan.json`:

- 64 public objects;
- 20 classes, 40 functions, and 4 values;
- 282 non-constructor methods and 180 documented attributes.

This contract maps that denominator into a smaller Rust-native set of stable
abstractions. A many-to-one mapping is allowed only when the shared Rust API
preserves every tested scientific behavior. A name appearing in a mapping
table is not evidence of implementation; completion is defined in `TESTS.md`.

## Public modules

```text
quspin
├── basis       built-in/general/composite/user bases, symmetry and projection
├── operator    parsed terms, stored and matrix-free operators, time dependence
├── solve       dense/partial eigensolvers, Lanczos, exponential action, evolve
├── dynamics    Floquet and spectral/dynamical response
├── measure     expectation, entropy, ensembles and time-series observables
├── block       block-diagonal construction and delayed block operations
├── archive     safe versioned operator persistence
├── workflow    subspace tracking and Lindblad generators
└── error       stable recoverable error categories
```

`basis::Basis` and `operator::LinearOperator` are the two required narrow
interfaces. Additional traits must represent genuine mathematical capability,
not storage or model names.

## Core state and basis boundary

The existing `Basis` contract is preserved:

```rust
pub trait Basis: Send + Sync {
    type State: Copy + Eq + Send + Sync;

    fn len(&self) -> usize;
    fn state(&self, index: usize) -> Result<Self::State>;
    fn index(&self, state: Self::State) -> Result<usize>;
    fn apply_local(
        &self,
        state: Self::State,
        operator: &str,
        sites: &[usize],
    ) -> Result<Option<(Self::State, Complex64)>>;
}
```

The full package adds these capability objects:

```rust
pub trait SymmetryMap<State>: Send + Sync {
    fn period(&self) -> usize;
    fn apply(&self, state: State) -> Result<(State, Complex64)>;
}

pub struct SymmetryReducer<State> { /* finite map, character and block order */ }

pub struct BasisProjector { /* source and reduced shapes plus sparse action */ }
impl LinearOperator for BasisProjector { /* projected action */ }

pub enum StateStorage {
    U128,
    U256,
    U1024,
    U4096,
    U16384,
}
```

Required built-in types are:

```rust
SpinBasis1D             SpinBasisGeneral
BosonBasis1D            BosonBasisGeneral
SpinlessFermionBasis1D  SpinlessFermionBasisGeneral
SpinfulFermionBasis1D   SpinfulFermionBasisGeneral
PhotonBasis             TensorBasis
UserBasis<State>
```

Builders validate particle sectors, local dimensions, map periods, character
compatibility, and state-backend capacity before construction. General and
user bases support deferred construction and explicit materialization.

The basis module must expose physical operations equivalent to:

- state/index round trips and integer/occupation conversion;
- representative, normalization, and symmetry amplitude;
- projection to and lifting from reduced sectors;
- partial trace and entanglement entropy;
- local in-place action and bra/ket transition tables;
- cross-sector action between compatible source and target bases;
- Hermiticity, particle-conservation, and symmetry checks.

## Linear operators and time dependence

The existing stored/matrix-free boundary is preserved:

```rust
pub trait LinearOperator: Send + Sync {
    fn shape(&self) -> (usize, usize);
    fn format(&self) -> MatrixFormat;
    fn apply(&self, input: &[Complex64], output: &mut [Complex64]) -> Result<()>;
    fn shifted_solver(&self, shift: f64)
        -> Result<Option<Box<dyn ShiftedLinearSolver>>>;
}

pub trait TimeDependentOperator: Send + Sync {
    fn shape(&self) -> (usize, usize);
    fn apply_at(
        &self,
        time: f64,
        input: &[Complex64],
        output: &mut [Complex64],
    ) -> Result<()>;
}
```

The required operator objects are:

```rust
OperatorSpec
DynamicTerm
OperatorBuilder<Source, Target>
Operator
Hamiltonian<Kind = Static>
QuantumLinearOperator
QuantumOperator
ExpOp
```

`OperatorBuilder::on` constructs square maps and `OperatorBuilder::between`
constructs rectangular maps. All formats use the same parsed local actions:

```rust
pub enum MatrixFormat { Dense, Csc, Csr, Dia, MatrixFree }
```

`Hamiltonian<Static>` implements both `LinearOperator` and
`TimeDependentOperator`. `Hamiltonian<Dynamic>` implements
`TimeDependentOperator`; evaluating it at a selected time produces a static
operator when materialization is requested. This type-state distinction keeps
static solver calls honest without duplicating term semantics.

`QuantumOperator` evaluates a named parameter dictionary without changing
basis or term semantics. Missing parameters follow the documented Python
default behavior and are tested explicitly.

The module provides addition, subtraction, scalar multiplication, products,
integer powers, commutator, anticommutator, expectation value, matrix element,
variance, projection, rotation, dense/sparse conversion, and diagonal access
where mathematically defined. Shape and time-dependence errors are explicit.

## Solver boundary

```rust
pub fn eigh(op: &impl LinearOperator, options: EighOptions)
    -> Result<Eigensystem>;

pub fn eigsh(op: &impl LinearOperator, options: EigshOptions)
    -> Result<Eigensystem>;

pub fn lanczos_full(op: &impl LinearOperator, initial: &[Complex64],
    options: LanczosOptions) -> Result<LanczosDecomposition>;

pub fn lanczos_iter<'a>(op: &'a impl LinearOperator,
    initial: &'a [Complex64], options: LanczosOptions)
    -> Result<impl Iterator<Item = Result<LanczosVector>> + 'a>;

pub fn expm_multiply(op: &impl LinearOperator, initial: &[Complex64],
    options: EvolutionOptions) -> Result<StateTrajectory>;

pub fn evolve<G: EvolutionGenerator>(generator: &G,
    initial: StateBatchRef<'_>, options: EvolutionOptions)
    -> Result<StateTrajectory>;
```

`StateBatchRef` represents either one state vector or a column-major batch.
`EvolutionGenerator` is implemented by static linear operators, explicitly
time-dependent operators, and callable right-hand sides with a documented
complex or stacked-real state convention.

The full `eigsh` target set covers algebraic ends, magnitude ends, both ends,
and real shift targeting. It accepts optional initial vectors, Krylov controls,
iteration caps, tolerances, deterministic seeds, and backend-specific options
through a typed extension that cannot silently change core semantics.

Dense and iterative Hermitian solvers support real-symmetric and complex-
Hermitian operators. Results include eigenvalues, eigenvectors when requested,
residuals, iteration counts, and convergence status.

Evolution accepts a scalar or ordered time grid, a state or batch of states,
Schrodinger and Liouville-von Neumann modes, time-independent and explicitly
time-dependent generators, and iterator or collected output. Norm or trace
preservation is a verified property when the generator implies it.

## Dynamics, block, measurement, and workflow APIs

Required dynamics:

```rust
Floquet
FloquetTimeVector
spectral_function
dynamical_correlator
```

`Floquet` supports piecewise and callable drives, stroboscopic action, optional
full-unitary construction, quasienergies/eigenvectors, and effective
Hamiltonians on small systems. `FloquetTimeVector` defines period, cycle, and
within-cycle coordinates without accumulating floating-point drift.

Required block tools:

```rust
BlockOps
block_diag_hamiltonian
```

Blocks may be computed eagerly or on demand. Static and dynamic block results
must reconstruct the corresponding full-space result within tolerance.

Required measurements and helpers:

```rust
expectation
matrix_element
quantum_fluctuation
partial_trace
entanglement_entropy
observables_vs_time
ed_state_vs_time
diagonal_ensemble
kl_divergence
mean_level_spacing
project_operator
array_to_states
states_to_array
```

Required thermal and Lanczos helpers:

```rust
ftlm_static_iteration
ltlm_static_iteration
expm_lanczos
linear_combination_qt
```

Existing workflow extensions remain public:

```rust
Subspace
subspace_fidelity
track_eigenspaces
LindbladGenerator
```

These extensions do not substitute for any frozen Python object; they are
additional workflow-derived capabilities.

## Archive contract

```rust
pub fn save_zip(path: impl AsRef<Path>, entries: &OperatorArchive)
    -> Result<()>;

pub fn load_zip(path: impl AsRef<Path>) -> Result<OperatorArchive>;
```

The archive preserves entry names, shapes, element types, dense or sparse
values, and sparse index structure. Round trips are lossless within the stored
numeric type. Python interoperability is verified in both directions. Loading
does not execute serialized code or callbacks.

## Error contract

All recoverable public failures return `Result<T, QmbedError>`. Stable
categories cover invalid operators/couplings/sites, unavailable states or
sectors, incompatible symmetry, non-Hermitian input, invalid options,
dimension mismatch, rank deficiency, unsupported backends/formats, archive
errors, and solver non-convergence with structured metadata.

Panics are reserved for violated internal invariants. Allocation failure and
process termination are outside the recoverable contract.

## Frozen Python object mapping

The status words below describe the required work from the current Rust crate:
`extend` means a partial analogue exists; `add` means no direct public analogue
exists; `preserve` means the current workflow API already supplies the planned
boundary but still needs full-package tests.

### Basis namespace — 28 objects

| Python object | Planned Rust boundary | Work |
|---|---|---|
| `basis_int_to_python_int` | `basis::state_to_biguint` | add |
| `basis_ones` | `basis::StateBits::ones` | add |
| `basis_zeros` | `basis::StateBits::zeros` | add |
| `bitwise_and` | `StateBits::bitand` and compatibility helper | add |
| `bitwise_leftshift` | `StateBits::shl` and compatibility helper | add |
| `bitwise_not` | `StateBits::not` and compatibility helper | add |
| `bitwise_or` | `StateBits::bitor` and compatibility helper | add |
| `bitwise_rightshift` | `StateBits::shr` and compatibility helper | add |
| `bitwise_xor` | `StateBits::bitxor` and compatibility helper | add |
| `python_int_to_basis_int` | `basis::state_from_biguint` | add |
| `uint256` | `basis::U256` | add |
| `uint1024` | `basis::U1024` | add |
| `uint4096` | `basis::U4096` | add |
| `uint16384` | `basis::U16384` | add |
| `get_basis_type` | `StateStorage::for_sites` | add |
| `spin_basis_1d` | `SpinBasis1D` | extend |
| `boson_basis_1d` | `BosonBasis1D` | extend |
| `spinless_fermion_basis_1d` | `SpinlessFermionBasis1D` | extend |
| `spinful_fermion_basis_1d` | `SpinfulFermionBasis1D` | extend |
| `spin_basis_general` | `SpinBasisGeneral` | add |
| `boson_basis_general` | `BosonBasisGeneral` | add |
| `spinless_fermion_basis_general` | `SpinlessFermionBasisGeneral` | add |
| `spinful_fermion_basis_general` | `SpinfulFermionBasisGeneral` | add |
| `user_basis` | `UserBasis<State>` plus symmetry/deferred capabilities | extend |
| `tensor_basis` | `TensorBasis` | add |
| `photon_basis` | `PhotonBasis` | add |
| `photon_Hspace_dim` | `PhotonBasis::dimension` | add |
| `coherent_state` | `basis::coherent_state` | add |

### Operator namespace — 12 objects

| Python object | Planned Rust boundary | Work |
|---|---|---|
| `hamiltonian` | `Hamiltonian` over `Operator` and `DynamicTerm` | extend |
| `quantum_LinearOperator` | `QuantumLinearOperator` / `LinearOperator` | extend |
| `quantum_operator` | `QuantumOperator` | add |
| `exp_op` | `ExpOp` | add |
| `commutator`, `anti_commutator` | `operator::{commutator,anticommutator}` | add |
| `isexp_op` | `OperatorKind` and compatibility helper | add |
| `ishamiltonian` | `OperatorKind` and compatibility helper | add |
| `isquantum_LinearOperator` | `OperatorKind` and compatibility helper | add |
| `isquantum_operator` | `OperatorKind` and compatibility helper | add |
| `save_zip`, `load_zip` | `archive::{save_zip,load_zip}` | add |

### Tools namespace — 24 objects

| Python object | Planned Rust boundary | Work |
|---|---|---|
| `Floquet` | `dynamics::Floquet` | extend |
| `Floquet_t_vec` | `dynamics::FloquetTimeVector` | add |
| `evolve` | `solve::evolve` | extend |
| `ExpmMultiplyParallel`, `expm_multiply_parallel` | `solve::ExpmPlan`, `expm_multiply` | add |
| `ED_state_vs_time` | `measure::ed_state_vs_time` | add |
| `block_ops`, `block_diag_hamiltonian` | `block::{BlockOps,block_diag_hamiltonian}` | add |
| `lanczos_full`, `lanczos_iter` | `solve::{lanczos_full,lanczos_iter}` | add |
| `expm_lanczos`, `lin_comb_Q_T` | `solve::{expm_lanczos,linear_combination_qt}` | add |
| `FTLM_static_iteration`, `LTLM_static_iteration` | thermal Lanczos helpers | add |
| `diag_ensemble`, `ent_entropy`, `obs_vs_time` | `measure` functions | add |
| `KL_div`, `mean_level_spacing` | `measure` statistics | add |
| `array_to_ints`, `ints_to_array` | `basis` state-array conversion | add |
| `get_matvec_function`, `matvec` | `LinearOperator` adapter helpers | add |
| `project_op` | `measure::project_operator` | add |

## Compatibility rules

- Current zero-based Rust site indexing remains stable.
- Current `OperatorBuilder::on/between`, `MatrixFormat`, and `LinearOperator`
  remain source-compatible unless a documented major-version migration is
  unavoidable.
- Enumeration order is not a physical contract; deterministic ordering is
  required for reproducibility, while tests compare through state lookup or
  projectors.
- Python-compatible helpers may live under `compat`, but core scientific
  operations must not require the compatibility module.
- Unsupported valid inputs are incomplete capability, not `InvalidOptions`.
  They may return `UnsupportedBackend` during staged development but cannot do
  so in a row marked complete.

## Completion gate

The full-package contract is complete only when:

1. every mapping row has a public Rust target and a public behavior test;
2. every required Python member is assigned to one Rust semantic test family;
3. independent Python-oracle and Python-independent property tests pass;
4. all twelve paper workflows pass through only the public candidate crate;
5. scale tests show that sparse, symmetry, cross-sector, and matrix-free paths
   avoid forbidden dense intermediates;
6. documentation and machine-readable denominators cannot shrink silently.
