//! Thin verification adapter from the frozen contract to QMBED.

use qmbed::basis::{
    Basis as PublicBasis, BosonBasis1D, SpinBasis1D, SpinfulFermionBasis1D, SpinlessFermionBasis1D,
    UserBasis,
};
use qmbed::dynamics::{spectral_function, SpectrumOptions as PublicSpectrumOptions};
use qmbed::measure::{subspace_fidelity, Subspace};
use qmbed::operator::{
    AssemblyChecks as PublicAssemblyChecks, Coupling as PublicCoupling,
    LinearOperator as PublicLinearOperator, MatrixFormat as PublicMatrixFormat, Operator,
    OperatorBuilder, OperatorTerm as PublicOperatorTerm,
};
use qmbed::solve::{
    eigsh, evolve, Eigensystem, EigshOptions as PublicEigshOptions,
    EvolutionOptions as PublicEvolutionOptions, SpectrumTarget as PublicSpectrumTarget,
};
use qmbed::{Complex64, QmbedError};

use crate::api::{
    BasisHandle, BasisSpec, EigshOptions, EvolutionOptions, HamiltonianOptions, LinearOperator,
    MatrixFormat, OperatorTerm, QmbedApi, SpectrumOptions, SpectrumTarget,
};

pub enum CandidateBasis {
    Spin(SpinBasis1D),
    Boson(BosonBasis1D),
    SpinlessFermion(SpinlessFermionBasis1D),
    SpinfulFermion(SpinfulFermionBasis1D),
    User(UserBasis<u128>),
}

impl CandidateBasis {
    fn as_basis(&self) -> &dyn PublicBasis<State = u128> {
        match self {
            Self::Spin(basis) => basis,
            Self::Boson(basis) => basis,
            Self::SpinlessFermion(basis) => basis,
            Self::SpinfulFermion(basis) => basis,
            Self::User(basis) => basis,
        }
    }
}

impl BasisHandle for CandidateBasis {
    type Error = QmbedError;

    fn dimension(&self) -> usize {
        self.as_basis().len()
    }

    fn state_at(&self, index: usize) -> Result<u128, Self::Error> {
        self.as_basis().state(index)
    }

    fn state_index(&self, state: u128) -> Result<usize, Self::Error> {
        self.as_basis().index(state)
    }
}

pub struct CandidateOperator(Operator);

impl LinearOperator for CandidateOperator {
    type Scalar = Complex64;
    type Error = QmbedError;

    fn shape(&self) -> (usize, usize) {
        self.0.shape()
    }

    fn format(&self) -> MatrixFormat {
        match self.0.format() {
            PublicMatrixFormat::Dense => MatrixFormat::Dense,
            PublicMatrixFormat::Csc => MatrixFormat::Csc,
            PublicMatrixFormat::Csr => MatrixFormat::Csr,
            PublicMatrixFormat::Dia => MatrixFormat::Dia,
            PublicMatrixFormat::MatrixFree => MatrixFormat::MatrixFree,
        }
    }

    fn apply(
        &self,
        input: &[Self::Scalar],
        output: &mut [Self::Scalar],
    ) -> Result<(), Self::Error> {
        self.0.apply(input, output)
    }
}

#[derive(Default)]
pub struct QmbedAdapter;

fn checked_count(value: i32, label: &str) -> Result<usize, QmbedError> {
    usize::try_from(value)
        .map_err(|_| QmbedError::InvalidSector(format!("{label} must be nonnegative")))
}

fn public_format(format: MatrixFormat) -> PublicMatrixFormat {
    match format {
        MatrixFormat::Dense => PublicMatrixFormat::Dense,
        MatrixFormat::Csc => PublicMatrixFormat::Csc,
        MatrixFormat::Csr => PublicMatrixFormat::Csr,
        MatrixFormat::Dia => PublicMatrixFormat::Dia,
        MatrixFormat::MatrixFree => PublicMatrixFormat::MatrixFree,
    }
}

fn public_terms(terms: &[OperatorTerm]) -> Result<Vec<PublicOperatorTerm>, QmbedError> {
    terms
        .iter()
        .map(|term| {
            PublicOperatorTerm::new(
                &term.operator,
                term.couplings.iter().map(|coupling| PublicCoupling {
                    coefficient: Complex64::new(coupling.coefficient.re, coupling.coefficient.im),
                    sites: coupling.sites.clone(),
                }),
            )
        })
        .collect()
}

fn build_operator<B>(
    basis: &B,
    terms: &[PublicOperatorTerm],
    options: HamiltonianOptions,
) -> Result<Operator, QmbedError>
where
    B: PublicBasis<State = u128>,
{
    OperatorBuilder::on(basis)
        .terms(terms.iter().cloned())
        .checks(PublicAssemblyChecks {
            hermiticity: options.checks.hermiticity,
            particle_conservation: options.checks.particle_conservation,
            symmetry_compatibility: options.checks.symmetry_compatibility,
        })
        .build(public_format(options.format))
}

fn user_basis(
    sites: usize,
    states: &[u128],
    allowed_operators: &[String],
) -> Result<UserBasis<u128>, QmbedError> {
    let mut builder = UserBasis::builder(sites).states(states.iter().copied());
    for operator in allowed_operators {
        builder = match operator.as_str() {
            "I" => builder.operator('I', |state, _| Ok(Some((state, Complex64::new(1.0, 0.0))))),
            "x" => builder.operator('x', |state, site| {
                Ok(Some((state ^ (1_u128 << site), Complex64::new(1.0, 0.0))))
            }),
            "z" => builder.operator('z', |state, site| {
                let value = if state & (1_u128 << site) == 0 {
                    -1.0
                } else {
                    1.0
                };
                Ok(Some((state, Complex64::new(value, 0.0))))
            }),
            name => return Err(QmbedError::InvalidOperator(name.to_string())),
        };
    }
    builder.build()
}

impl QmbedApi for QmbedAdapter {
    type Scalar = Complex64;
    type Error = QmbedError;
    type Basis = CandidateBasis;
    type Operator = CandidateOperator;
    type State = Vec<Complex64>;
    type Eigensystem = Eigensystem;

    fn basis(&self, spec: &BasisSpec) -> Result<Self::Basis, Self::Error> {
        match spec {
            BasisSpec::Spin {
                sites,
                spin_twice,
                magnetization,
                momentum,
                parity,
                pauli,
            } => {
                let mut builder = SpinBasis1D::builder(*sites)
                    .spin_twice(*spin_twice)
                    .pauli(*pauli);
                if let Some(count) = magnetization {
                    builder = builder.magnetization(checked_count(*count, "magnetization")?);
                }
                if let Some(sector) = momentum {
                    builder = builder.momentum(*sector);
                }
                if let Some(sector) = parity {
                    builder = builder.parity(*sector);
                }
                builder.build().map(CandidateBasis::Spin)
            }
            BasisSpec::Boson {
                sites,
                particles,
                states_per_site,
            } => {
                let mut builder = BosonBasis1D::builder(*sites, *states_per_site);
                if let Some(count) = particles {
                    builder = builder.particles(*count);
                }
                builder.build().map(CandidateBasis::Boson)
            }
            BasisSpec::SpinlessFermion {
                sites,
                particles,
                momentum,
            } => {
                let mut builder = SpinlessFermionBasis1D::builder(*sites);
                if let Some(count) = particles {
                    builder = builder.particles(*count);
                }
                if let Some(sector) = momentum {
                    builder = builder.momentum(*sector);
                }
                builder.build().map(CandidateBasis::SpinlessFermion)
            }
            BasisSpec::SpinfulFermion {
                sites,
                particles_up,
                particles_down,
            } => {
                let mut builder = SpinfulFermionBasis1D::builder(*sites);
                if let Some(count) = particles_up {
                    builder = builder.particles_up(*count);
                }
                if let Some(count) = particles_down {
                    builder = builder.particles_down(*count);
                }
                builder.build().map(CandidateBasis::SpinfulFermion)
            }
            BasisSpec::User {
                sites,
                states,
                allowed_operators,
            } => user_basis(*sites, states, allowed_operators).map(CandidateBasis::User),
            BasisSpec::Photon { .. } | BasisSpec::Tensor { .. } => Err(
                QmbedError::UnsupportedBackend("basis is outside the frozen Rust surface".into()),
            ),
        }
    }

    fn hamiltonian(
        &self,
        basis: &Self::Basis,
        terms: &[OperatorTerm],
        options: HamiltonianOptions,
    ) -> Result<Self::Operator, Self::Error> {
        let terms = public_terms(terms)?;
        let operator = match basis {
            CandidateBasis::Spin(basis) => build_operator(basis, &terms, options),
            CandidateBasis::Boson(basis) => build_operator(basis, &terms, options),
            CandidateBasis::SpinlessFermion(basis) => build_operator(basis, &terms, options),
            CandidateBasis::SpinfulFermion(basis) => build_operator(basis, &terms, options),
            CandidateBasis::User(basis) => build_operator(basis, &terms, options),
        }?;
        Ok(CandidateOperator(operator))
    }

    fn eigsh(
        &self,
        operator: &Self::Operator,
        _initial: Option<&Self::State>,
        options: &EigshOptions,
    ) -> Result<Self::Eigensystem, Self::Error> {
        let target = match options.target {
            SpectrumTarget::SmallestAlgebraic => PublicSpectrumTarget::SmallestAlgebraic,
            SpectrumTarget::LargestAlgebraic => PublicSpectrumTarget::LargestAlgebraic,
            SpectrumTarget::LargestMagnitude => PublicSpectrumTarget::LargestMagnitude,
            SpectrumTarget::Shift(value) => PublicSpectrumTarget::Shift(value),
        };
        eigsh(
            &operator.0,
            PublicEigshOptions {
                eigenpairs: options.eigenpairs,
                target,
                krylov_dimension: options.krylov_dimension,
                tolerance: options.tolerance,
                max_iterations: options.max_iterations,
                seed: 0,
            },
        )
    }

    fn evolve(
        &self,
        operator: &Self::Operator,
        initial: &Self::State,
        options: &EvolutionOptions,
    ) -> Result<Vec<Self::State>, Self::Error> {
        evolve(
            &operator.0,
            initial,
            PublicEvolutionOptions {
                times: options.times.clone(),
                krylov_dimension: options.krylov_dimension,
                tolerance: options.tolerance,
                max_substeps: 10_000,
                hamiltonian: true,
            },
        )
        .map(|trajectory| trajectory.states)
    }

    fn spectral_function(
        &self,
        operator: &Self::Operator,
        source: &Self::State,
        probe: &Self::Operator,
        options: &SpectrumOptions,
    ) -> Result<Vec<f64>, Self::Error> {
        spectral_function(
            &operator.0,
            source,
            &probe.0,
            PublicSpectrumOptions {
                frequencies: options.frequencies.clone(),
                reference_energy: options.reference_energy,
                broadening: options.broadening,
                krylov_dimension: options.krylov_dimension,
                tolerance: 1.0e-10,
            },
        )
    }

    fn subspace_fidelity(
        &self,
        left: &[Self::State],
        right: &[Self::State],
    ) -> Result<f64, Self::Error> {
        let ambient = left.first().or_else(|| right.first()).map_or(0, Vec::len);
        if left
            .iter()
            .chain(right)
            .any(|vector| vector.len() != ambient)
        {
            return Err(QmbedError::DimensionMismatch(
                "subspace vectors must share an ambient dimension".into(),
            ));
        }
        let left = Subspace::from_columns(
            ambient,
            left.len(),
            left.iter()
                .flat_map(|vector| vector.iter().copied())
                .collect(),
        )?;
        let right = Subspace::from_columns(
            ambient,
            right.len(),
            right
                .iter()
                .flat_map(|vector| vector.iter().copied())
                .collect(),
        )?;
        subspace_fidelity(&left, &right)
    }

    fn lindblad_evolve(
        &self,
        generator: &Self::Operator,
        initial_density: &Self::State,
        options: &EvolutionOptions,
    ) -> Result<Vec<Self::State>, Self::Error> {
        evolve(
            &generator.0,
            initial_density,
            PublicEvolutionOptions {
                times: options.times.clone(),
                krylov_dimension: options.krylov_dimension,
                tolerance: options.tolerance,
                max_substeps: 10_000,
                hamiltonian: false,
            },
        )
        .map(|trajectory| trajectory.states)
    }
}
