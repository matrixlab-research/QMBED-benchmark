# quspin — Contract (frozen)

> Generated from frozen SpecBundle `sha256:a9093788a8b2687689113102daaef231e58d33f6188b74927e9b60c097ee1721` — do not hand-edit.
> Everything in this file is immutable for the duration of the task.

**Deliverable.** A standard, publishable Rust crate named
`quspin` using the 2024 edition. Keep the ordinary Cargo layout (`Cargo.toml`,
`src/lib.rs`, and modules below it). The table is the **required public
surface**: expose every listed item at its fully qualified path with compatible
signatures and behavior. Additional public conveniences are allowed but are
not judged. Internal representation, helper items, dependencies, ownership
details, and algorithm choice are yours.

## Comparison semantics

The `comparison` column tells you how each function's output is checked, so you
know the required accuracy:

- `num:atol:rtol` — numeric agreement `|got − expected| ≤ atol + rtol·|expected|`
  (the stated scalar or element type); `±Inf` must match exactly, `NaN` maps to
  `NaN`; elementwise for same-shape numeric arrays.
- `isequal` — exact structural or scalar equality.
- `str` — the printed representation matches.
- `custom:<name>` — a named comparison described where it is used.

Returning `Err` or panicking on a valid input is a failure. Panics are reserved for internal invariant violations.

## Symbols

| symbol | signatures | comparison | domain |
|---|---|---|---|
| `error::QuSpinError` | `pub enum QuSpinError { InvalidOperator, InvalidCoupling, InvalidSite, StateNotInBasis, InvalidSector, IncompatibleSymmetry, NonHermitian, InvalidOptions, DimensionMismatch, RankDeficient, UnsupportedBackend, NonConvergence { iterations: usize, residual: f64 } }`; `pub type Result<T> = std::result::Result<T, QuSpinError>` | `custom:error-category` | all invalid user inputs and recoverable solver or backend failures |
| `basis::Basis` | `pub trait Basis: Send + Sync { type State: Copy + Eq + Send + Sync; fn len(&self) -> usize; fn state(&self, index: usize) -> Result<Self::State>; fn index(&self, state: Self::State) -> Result<usize>; fn apply_local(&self, state: Self::State, operator: &str, sites: &[usize]) -> Result<Option<(Self::State, Complex64)>>; }` | `custom:basis-observation` | finite discrete Hilbert spaces whose states have a stable equality relation; state ordering is implementation-defined |
| `basis::SpinBasis1D` | `pub fn builder(sites: usize) -> SpinBasisBuilder` | `custom:basis-observation` | sites fit the selected state representation; spin_twice is positive; magnetization and symmetry labels describe a nonempty compatible sector |
| `basis::BosonBasis1D` | `pub fn builder(sites: usize, states_per_site: usize) -> BosonBasisBuilder` | `custom:basis-observation` | sites and states_per_site are positive; optional total particle count is representable under the local cutoff |
| `basis::SpinlessFermionBasis1D` | `pub fn builder(sites: usize) -> SpinlessFermionBasisBuilder` | `custom:basis-observation` | finite one-dimensional lattices with optional fixed particle number and compatible lattice symmetry |
| `basis::SpinfulFermionBasis1D` | `pub fn builder(sites: usize) -> SpinfulFermionBasisBuilder` | `custom:basis-observation` | finite lattices with independently optional up and down particle counts |
| `basis::UserBasis` | `pub fn builder<State>(sites: usize) -> UserBasisBuilder<State> where State: Copy + Eq + std::hash::Hash + Send + Sync` | `custom:basis-observation` | finite user-defined state spaces with deterministic filtering and registered local operator callbacks |
| `operator::OperatorTerm` | `pub fn new(operator: impl AsRef<str>, couplings: impl IntoIterator<Item = Coupling>) -> Result<Self>` | `custom:parsed-term` | operator strings registered by the chosen basis and couplings whose site arity equals the operator-string arity |
| `operator::Coupling` | `pub struct Coupling { pub coefficient: Complex64, pub sites: Vec<usize> }` | `isequal` | finite complex coefficients and zero-based site lists whose arity is validated by OperatorTerm |
| `operator::MatrixFormat` | `pub enum MatrixFormat { Dense, Csc, Csr, Dia, MatrixFree }` | `isequal` | explicit backend requests supported by the selected operator and algorithm |
| `operator::LinearOperator` | `pub trait LinearOperator: Send + Sync { fn shape(&self) -> (usize, usize); fn format(&self) -> MatrixFormat; fn apply(&self, input: &[Complex64], output: &mut [Complex64]) -> Result<()>; }` | `num:1e-12:1e-10` | finite complex linear maps, including rectangular maps between distinct basis sectors |
| `operator::OperatorBuilder` | `pub fn on<B: Basis>(basis: &B) -> OperatorBuilder<'_, B, B>`; `pub fn between<S, T>(source: &S, target: &T) -> OperatorBuilder<'_, S, T> where S: Basis, T: Basis<State = S::State>`; `pub fn build(self, format: MatrixFormat) -> Result<Operator>` | `num:1e-12:1e-10` | valid parsed terms whose sites lie in the physical lattice and whose local actions map from the source basis to the target basis |
| `solve::EigshOptions` | `pub struct EigshOptions { pub eigenpairs: usize, pub target: SpectrumTarget, pub krylov_dimension: Option<usize>, pub tolerance: f64, pub max_iterations: usize, pub seed: u64 }`; `pub enum SpectrumTarget { SmallestAlgebraic, LargestAlgebraic, LargestMagnitude, Shift(f64) }` | `isequal` | positive eigenpair count, tolerance, and maximum iterations; optional Krylov dimension exceeds eigenpair count; shift is finite |
| `solve::eigsh` | `pub fn eigsh(op: &impl LinearOperator, options: EigshOptions) -> Result<Eigensystem>` | `num:1e-09:1e-08` | square Hermitian complex operators; requested eigenpair count is positive and smaller than the dimension; target is extremal, magnitude, or real shift |
| `solve::EvolutionOptions` | `pub struct EvolutionOptions { pub times: Vec<f64>, pub krylov_dimension: usize, pub tolerance: f64, pub max_substeps: usize, pub hamiltonian: bool }` | `isequal` | nonempty finite nondecreasing time grid, positive Krylov dimension, tolerance, and substep cap |
| `solve::evolve` | `pub fn evolve(op: &impl LinearOperator, initial: &[Complex64], options: EvolutionOptions) -> Result<StateTrajectory>` | `num:1e-09:1e-08` | square operators, an initial state of matching dimension, and a finite scalar time or finite ordered time grid |
| `dynamics::DriveStep` | `pub struct DriveStep { pub hamiltonian: std::sync::Arc<dyn LinearOperator>, pub duration: f64 }` | `custom:drive-step` | square Hermitian operators and finite nonnegative durations |
| `dynamics::Floquet` | `pub fn new(steps: impl IntoIterator<Item = DriveStep>) -> Result<Self>`; `pub fn apply_period(&self, input: &[Complex64], output: &mut [Complex64]) -> Result<()>` | `num:1e-09:1e-08` | nonempty ordered drive steps with square Hermitian Hamiltonians of equal dimension and finite nonnegative durations |
| `dynamics::SpectrumOptions` | `pub struct SpectrumOptions { pub frequencies: Vec<f64>, pub reference_energy: f64, pub broadening: f64, pub krylov_dimension: usize, pub tolerance: f64 }` | `isequal` | nonempty finite frequency grid, finite reference energy, positive finite broadening, Krylov dimension, and tolerance |
| `dynamics::spectral_function` | `pub fn spectral_function(h_target: &impl LinearOperator, source: &[Complex64], probe: &impl LinearOperator, options: SpectrumOptions) -> Result<Vec<f64>>` | `num:1e-08:1e-07` | a square target-sector Hamiltonian, a source-sector state, a probe whose columns match the source and rows match the target, positive broadening, and finite frequencies |
| `measure::Subspace` | `pub fn from_columns(ambient_dimension: usize, rank: usize, column_major_vectors: Vec<Complex64>) -> Result<Self>` | `num:1e-12:1e-10` | nonempty finite column sets with ambient_dimension times rank entries and full numerical column rank |
| `measure::subspace_fidelity` | `pub fn subspace_fidelity(left: &Subspace, right: &Subspace) -> Result<f64>` | `num:1e-12:1e-10` | two nonempty finite-dimensional subspaces represented by equal-ambient-dimension column sets; columns need not already be orthonormal |
| `workflow::LindbladGenerator` | `pub fn new(hamiltonian: std::sync::Arc<dyn LinearOperator>, jumps: Vec<std::sync::Arc<dyn LinearOperator>>) -> Result<Self>` | `num:1e-09:1e-08` | square Hamiltonian and jump operators of a common Hilbert dimension; density matrices use documented column-major vectorization |

## How you are judged

- Your package is tested by a **private verification suite** (held-out golden
  values in the edge/overflow/underflow regions, property checks over
  randomized inputs, and real multi-call integration usage). You never see that
  suite; the visible examples in TESTS.md are a floor, not the gate.
- Feedback is the **pass/fail of the verification pipeline**, mirrored onto your
  merge request. The merge request can merge only when it is green.
- Correctness rests on the spec: derive behavior from MOTIVATION.md and the
  comparison tolerances, not from guessing at the hidden cases. Passing the
  visible examples does not imply passing the gate.
- After 5 failed verification runs on a branch, the attempt is
  regenerated from the spec rather than patched further — so aim for a coherent
  implementation from the definitions, not trial-and-error against the gate.
- Scope: every listed Rust item is required and judged; additional public conveniences are permitted but are outside the gate.
