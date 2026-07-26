//! Candidate-facing API adapter contract.
//!
//! The concrete Rust package may use generic basis and matrix types internally.
//! The public verifier only requires that an adapter can express this narrow
//! set of operations without exposing implementation details.

use std::error::Error;

/// Sparse/dense materialization requested by a caller.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatrixFormat {
    Dense,
    Csc,
    Csr,
    Dia,
    MatrixFree,
}

impl MatrixFormat {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Dense => "dense",
            Self::Csc => "csc",
            Self::Csr => "csr",
            Self::Dia => "dia",
            Self::MatrixFree => "matrix_free",
        }
    }
}

/// Typed basis construction requests used by verification adapters.
#[derive(Clone, Debug, PartialEq)]
pub enum BasisSpec {
    Spin {
        sites: usize,
        spin_twice: u16,
        magnetization: Option<i32>,
        momentum: Option<i32>,
        parity: Option<i8>,
        pauli: bool,
    },
    Boson {
        sites: usize,
        particles: Option<usize>,
        states_per_site: usize,
    },
    SpinlessFermion {
        sites: usize,
        particles: Option<usize>,
        momentum: Option<i32>,
    },
    SpinfulFermion {
        sites: usize,
        particles_up: Option<usize>,
        particles_down: Option<usize>,
    },
    Photon {
        cutoff: usize,
    },
    Tensor {
        factors: Vec<BasisSpec>,
    },
    User {
        sites: usize,
        states: Vec<u128>,
        allowed_operators: Vec<String>,
    },
}

/// A complex coefficient without committing the verification harness to one
/// array or sparse-matrix crate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexCoefficient {
    pub re: f64,
    pub im: f64,
}

impl From<f64> for ComplexCoefficient {
    fn from(value: f64) -> Self {
        Self { re: value, im: 0.0 }
    }
}

/// One coefficient and its zero-based sites.
#[derive(Clone, Debug, PartialEq)]
pub struct Coupling {
    pub coefficient: ComplexCoefficient,
    pub sites: Vec<usize>,
}

/// Parsed once, applied many times while assembling or acting matrix-free.
#[derive(Clone, Debug, PartialEq)]
pub struct OperatorSpec {
    pub operator: String,
    pub couplings: Vec<Coupling>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AssemblyChecks {
    pub hermiticity: bool,
    pub particle_conservation: bool,
    pub symmetry_compatibility: bool,
}

impl Default for AssemblyChecks {
    fn default() -> Self {
        Self {
            hermiticity: true,
            particle_conservation: true,
            symmetry_compatibility: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct HamiltonianOptions {
    pub format: MatrixFormat,
    pub checks: AssemblyChecks,
}

impl Default for HamiltonianOptions {
    fn default() -> Self {
        Self {
            format: MatrixFormat::Csc,
            checks: AssemblyChecks::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpectrumTarget {
    SmallestAlgebraic,
    LargestAlgebraic,
    LargestMagnitude,
    Shift(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub struct EigshOptions {
    pub eigenpairs: usize,
    pub target: SpectrumTarget,
    pub krylov_dimension: Option<usize>,
    pub tolerance: f64,
    pub max_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EvolutionOptions {
    pub times: Vec<f64>,
    pub krylov_dimension: usize,
    pub tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumOptions {
    pub frequencies: Vec<f64>,
    pub reference_energy: f64,
    pub broadening: f64,
    pub krylov_dimension: usize,
}

/// Minimal basis behavior required by the verifier.
#[allow(clippy::missing_errors_doc)]
pub trait BasisHandle {
    type Error: Error + Send + Sync + 'static;

    fn dimension(&self) -> usize;
    fn state_at(&self, index: usize) -> Result<u128, Self::Error>;
    fn state_index(&self, state: u128) -> Result<usize, Self::Error>;
}

/// Rectangular-capable narrow waist shared by stored and matrix-free maps.
#[allow(clippy::missing_errors_doc)]
pub trait LinearOperator {
    type Scalar: Copy + Default + Send + Sync + 'static;
    type Error: Error + Send + Sync + 'static;

    fn shape(&self) -> (usize, usize);
    fn format(&self) -> MatrixFormat;
    fn apply(&self, input: &[Self::Scalar], output: &mut [Self::Scalar])
        -> Result<(), Self::Error>;
}
/// Verification-side adapter over the proposed public Rust API.
///
/// The future package is not required to implement this trait directly. The
/// verification repository can own a thin wrapper, avoiding an orphan-rule coupling
/// between the public crate and the verifier.
#[allow(clippy::missing_errors_doc)]
pub trait QmbedApi {
    type Scalar: Copy + Default + Send + Sync + 'static;
    type Error: Error + Send + Sync + 'static;
    type Basis: BasisHandle<Error = Self::Error>;
    type Operator: LinearOperator<Scalar = Self::Scalar, Error = Self::Error>;
    type State: AsRef<[Self::Scalar]>;
    type Eigensystem;

    fn basis(&self, spec: &BasisSpec) -> Result<Self::Basis, Self::Error>;

    fn hamiltonian(
        &self,
        basis: &Self::Basis,
        terms: &[OperatorSpec],
        options: HamiltonianOptions,
    ) -> Result<Self::Operator, Self::Error>;

    fn eigsh(
        &self,
        operator: &Self::Operator,
        initial: Option<&Self::State>,
        options: &EigshOptions,
    ) -> Result<Self::Eigensystem, Self::Error>;

    fn evolve(
        &self,
        operator: &Self::Operator,
        initial: &Self::State,
        options: &EvolutionOptions,
    ) -> Result<Vec<Self::State>, Self::Error>;

    fn spectral_function(
        &self,
        operator: &Self::Operator,
        source: &Self::State,
        probe: &Self::Operator,
        options: &SpectrumOptions,
    ) -> Result<Vec<f64>, Self::Error>;

    fn subspace_fidelity(
        &self,
        left: &[Self::State],
        right: &[Self::State],
    ) -> Result<f64, Self::Error>;

    fn lindblad_evolve(
        &self,
        generator: &Self::Operator,
        initial_density: &Self::State,
        options: &EvolutionOptions,
    ) -> Result<Vec<Self::State>, Self::Error>;
}
