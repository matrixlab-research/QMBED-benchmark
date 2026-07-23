//! End-to-end paper workflows executed against the pinned public candidate.

use std::collections::BTreeSet;
use std::sync::Arc;

use nalgebra::{linalg::Schur, DMatrix};
use quspin::basis::{
    Basis, BosonBasis1D, SpinBasis1D, SpinfulFermionBasis1D, SpinlessFermionBasis1D, UserBasis,
};
use quspin::dynamics::{spectral_function, DriveStep, Floquet, SpectrumOptions};
use quspin::measure::{subspace_fidelity, Subspace};
use quspin::operator::{
    Coupling, LinearOperator, MatrixFormat as PublicMatrixFormat, OperatorBuilder, OperatorTerm,
};
use quspin::solve::{eigsh, evolve, EigshOptions, EvolutionOptions, SpectrumTarget};
use quspin::{Complex64, QuSpinError};

use crate::candidate::QuSpinAdapter;
use crate::{MatrixFormat, Observation, WorkflowBackend, WorkflowCase};

fn c(value: f64) -> Complex64 {
    Complex64::new(value, 0.0)
}

fn observation() -> Observation {
    Observation::new(MatrixFormat::Csc)
}

fn require_dimension(actual: usize, expected: usize) -> Result<(), QuSpinError> {
    if actual == expected {
        Ok(())
    } else {
        Err(QuSpinError::InvalidSector(format!(
            "paper workflow expected dimension {expected}, got {actual}"
        )))
    }
}

fn maximum_residual(residuals: &[f64]) -> f64 {
    residuals.iter().copied().fold(0.0_f64, f64::max)
}

fn state_norm(state: &[Complex64]) -> f64 {
    state.iter().map(Complex64::norm_sqr).sum::<f64>().sqrt()
}

fn usize_as_f64(value: usize) -> Result<f64, QuSpinError> {
    u32::try_from(value).map(f64::from).map_err(|_| {
        QuSpinError::InvalidOptions(format!("workflow count {value} exceeds exact f64 range"))
    })
}

fn periodic_heisenberg_terms(sites: usize) -> Result<[OperatorTerm; 3], QuSpinError> {
    let mut zz = Vec::with_capacity(sites);
    let mut forward = Vec::with_capacity(sites);
    let mut backward = Vec::with_capacity(sites);
    for site in 0..sites {
        let next = (site + 1) % sites;
        zz.push(Coupling::new(1.0, vec![site, next]));
        forward.push(Coupling::new(0.5, vec![site, next]));
        backward.push(Coupling::new(0.5, vec![site, next]));
    }
    Ok([
        OperatorTerm::new("zz", zz)?,
        OperatorTerm::new("+-", forward)?,
        OperatorTerm::new("-+", backward)?,
    ])
}

fn periodic_blockade_states(sites: usize) -> Vec<u128> {
    fn extend(
        site: usize,
        sites: usize,
        first_occupied: bool,
        previous_occupied: bool,
        state: u128,
        output: &mut Vec<u128>,
    ) {
        if site == sites {
            if !(first_occupied && previous_occupied) {
                output.push(state);
            }
            return;
        }
        extend(site + 1, sites, first_occupied, false, state, output);
        if !previous_occupied {
            extend(
                site + 1,
                sites,
                first_occupied || site == 0,
                true,
                state | (1_u128 << site),
                output,
            );
        }
    }

    let mut states = Vec::new();
    extend(0, sites, false, false, 0, &mut states);
    states.sort_unstable();
    states
}

fn mbl_shift_invert() -> Result<Observation, QuSpinError> {
    let sites = 14;
    let basis = SpinBasis1D::builder(sites).up(7).build()?;
    require_dimension(basis.len(), 3_432)?;
    let fields = [
        2.13, -1.77, 0.31, 3.24, -2.63, 0.82, 1.46, -3.17, 2.71, -0.54, 1.09, -2.28, 0.67, 2.94,
    ];
    let mut terms = periodic_heisenberg_terms(sites)?.to_vec();
    terms.push(OperatorTerm::new(
        "z",
        fields
            .into_iter()
            .enumerate()
            .map(|(site, field)| Coupling::new(field, vec![site])),
    )?);
    let hamiltonian = OperatorBuilder::on(&basis)
        .terms(terms)
        .build(PublicMatrixFormat::Csc)?;
    let result = eigsh(
        &hamiltonian,
        EigshOptions {
            eigenpairs: 6,
            target: SpectrumTarget::Shift(0.0),
            krylov_dimension: Some(32),
            tolerance: 1.0e-9,
            max_iterations: 5_000,
            seed: 47,
        },
    )?;
    let residual = result
        .residuals
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    Ok(observation().metric("residual", residual))
}

fn xxz_lanczos_quench() -> Result<Observation, QuSpinError> {
    let sites = 16;
    let basis = SpinBasis1D::builder(sites).up(8).build()?;
    require_dimension(basis.len(), 12_870)?;
    let bonds = 0..(sites - 1);
    let hamiltonian = OperatorBuilder::on(&basis)
        .terms([
            OperatorTerm::new(
                "+-",
                bonds
                    .clone()
                    .map(|site| Coupling::new(0.5, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "-+",
                bonds
                    .clone()
                    .map(|site| Coupling::new(0.5, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "zz",
                bonds.map(|site| Coupling::new(0.8, vec![site, site + 1])),
            )?,
        ])
        .build(PublicMatrixFormat::Csc)?;
    let neel = (0..sites)
        .step_by(2)
        .fold(0_u128, |state, site| state | (1_u128 << site));
    let mut initial = vec![c(0.0); basis.len()];
    initial[basis.index(neel)?] = c(1.0);
    let trajectory = evolve(
        &hamiltonian,
        &initial,
        EvolutionOptions {
            times: vec![0.7],
            krylov_dimension: 80,
            tolerance: 1.0e-10,
            max_substeps: 100,
            hamiltonian: true,
        },
    )?;
    let norm_error = (state_norm(&trajectory.states[0]) - 1.0).abs();
    Ok(observation().metric("norm_error", norm_error))
}

fn floquet_heating() -> Result<Observation, QuSpinError> {
    let sites = 9;
    let basis = SpinBasis1D::builder(sites).pauli(true).build()?;
    require_dimension(basis.len(), 512)?;
    let zz = OperatorBuilder::on(&basis)
        .term(OperatorTerm::new(
            "zz",
            (0..sites).map(|site| Coupling::new(0.9, vec![site, (site + 1) % sites])),
        )?)
        .build(PublicMatrixFormat::Csc)?;
    let x = OperatorBuilder::on(&basis)
        .term(OperatorTerm::new(
            "x",
            (0..sites).map(|site| Coupling::new(0.73, vec![site])),
        )?)
        .build(PublicMatrixFormat::Csc)?;
    let floquet = Floquet::new([
        DriveStep::new(Arc::new(zz), 0.17)?,
        DriveStep::new(Arc::new(x), 0.23)?,
    ])?;
    let dimension = basis.len();
    let mut column_major = vec![c(0.0); dimension * dimension];
    let mut input = vec![c(0.0); dimension];
    let mut output = vec![c(0.0); dimension];
    for column in 0..dimension {
        input.fill(c(0.0));
        input[column] = c(1.0);
        floquet.apply_period(&input, &mut output)?;
        for row in 0..dimension {
            column_major[row + column * dimension] = output[row];
        }
    }
    let unitary = DMatrix::from_column_slice(dimension, dimension, &column_major);
    let gram = unitary.adjoint() * &unitary;
    let unitarity_error = (gram - DMatrix::<Complex64>::identity(dimension, dimension)).norm()
        / usize_as_f64(dimension)?;
    let (_, triangular) = Schur::new(unitary).unpack();
    let phase_modulus_error = (0..dimension)
        .map(|index| (triangular[(index, index)].norm() - 1.0).abs())
        .fold(0.0_f64, f64::max);
    Ok(observation()
        .metric("unitarity_error", unitarity_error)
        .metric("phase_modulus_error", phase_modulus_error))
}

fn spinful_hubbard() -> Result<Observation, QuSpinError> {
    let sites = 8;
    let basis = SpinfulFermionBasis1D::builder(sites)
        .particles(4, 4)
        .build()?;
    require_dimension(basis.len(), 4_900)?;
    let bonds = 0..(sites - 1);
    let hamiltonian = OperatorBuilder::on(&basis)
        .terms([
            OperatorTerm::new(
                "+-|",
                bonds
                    .clone()
                    .map(|site| Coupling::new(-1.0, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "-+|",
                bonds
                    .clone()
                    .map(|site| Coupling::new(1.0, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "|+-",
                bonds
                    .clone()
                    .map(|site| Coupling::new(-1.0, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "|-+",
                bonds.map(|site| Coupling::new(1.0, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "n|n",
                (0..sites).map(|site| Coupling::new(4.0, vec![site, site])),
            )?,
        ])
        .build(PublicMatrixFormat::Csc)?;
    let result = eigsh(
        &hamiltonian,
        EigshOptions {
            eigenpairs: 6,
            target: SpectrumTarget::SmallestAlgebraic,
            krylov_dimension: Some(160),
            tolerance: 1.0e-9,
            max_iterations: 192,
            seed: 37,
        },
    )?;
    Ok(observation().metric("residual", maximum_residual(&result.residuals)))
}

fn interacting_ssh() -> Result<Observation, QuSpinError> {
    let sites = 16;
    let basis = SpinlessFermionBasis1D::builder(sites)
        .particles(8)
        .build()?;
    require_dimension(basis.len(), 12_870)?;
    let hopping = |site: usize| if site % 2 == 0 { 0.6 } else { 1.0 };
    let bonds = 0..(sites - 1);
    let hamiltonian = OperatorBuilder::on(&basis)
        .terms([
            OperatorTerm::new(
                "+-",
                bonds
                    .clone()
                    .map(|site| Coupling::new(-hopping(site), vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "-+",
                bonds
                    .clone()
                    .map(|site| Coupling::new(hopping(site), vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "nn",
                bonds.map(|site| Coupling::new(2.0, vec![site, site + 1])),
            )?,
        ])
        .build(PublicMatrixFormat::Csc)?;
    let result = eigsh(
        &hamiltonian,
        EigshOptions {
            eigenpairs: 6,
            target: SpectrumTarget::SmallestAlgebraic,
            krylov_dimension: Some(160),
            tolerance: 1.0e-9,
            max_iterations: 192,
            seed: 41,
        },
    )?;
    Ok(observation().metric("residual", maximum_residual(&result.residuals)))
}

fn translation_sector_xxz() -> Result<Observation, QuSpinError> {
    let basis = SpinBasis1D::builder(18).up(9).momentum(0).build()?;
    require_dimension(basis.len(), 2_704)?;
    let hamiltonian = OperatorBuilder::on(&basis)
        .terms(periodic_heisenberg_terms(18)?)
        .build(PublicMatrixFormat::Csc)?;
    let result = eigsh(
        &hamiltonian,
        EigshOptions {
            eigenpairs: 4,
            target: SpectrumTarget::SmallestAlgebraic,
            krylov_dimension: Some(96),
            tolerance: 1.0e-8,
            max_iterations: 128,
            seed: 31,
        },
    )?;
    Ok(observation().metric("residual", maximum_residual(&result.residuals)))
}

fn tfim_fidelity_scan() -> Result<Observation, QuSpinError> {
    let sites = 16;
    let basis = SpinBasis1D::builder(sites).pauli(true).build()?;
    require_dimension(basis.len(), 65_536)?;
    let mut subspaces = Vec::new();
    let mut residual = 0.0_f64;
    for (field_index, field) in [0.8, 0.9, 1.0, 1.1, 1.2].into_iter().enumerate() {
        let hamiltonian = OperatorBuilder::on(&basis)
            .terms([
                OperatorTerm::new(
                    "zz",
                    (0..sites).map(|site| Coupling::new(-1.0, vec![site, (site + 1) % sites])),
                )?,
                OperatorTerm::new(
                    "x",
                    (0..sites).map(|site| Coupling::new(-field, vec![site])),
                )?,
            ])
            .build(PublicMatrixFormat::Csc)?;
        let result = eigsh(
            &hamiltonian,
            EigshOptions {
                eigenpairs: 2,
                target: SpectrumTarget::SmallestAlgebraic,
                krylov_dimension: Some(100),
                tolerance: 1.0e-9,
                max_iterations: 128,
                seed: 43 + field_index as u64,
            },
        )?;
        residual = residual.max(maximum_residual(&result.residuals));
        subspaces.push(Subspace::from_columns(
            basis.len(),
            2,
            result.eigenvectors.into_iter().flatten().collect(),
        )?);
    }
    let fidelities = subspaces
        .windows(2)
        .map(|pair| subspace_fidelity(&pair[0], &pair[1]))
        .collect::<Result<Vec<_>, _>>()?;
    let (minimum_index, minimum_fidelity) = fidelities
        .iter()
        .enumerate()
        .min_by(|left, right| left.1.total_cmp(right.1))
        .ok_or_else(|| QuSpinError::InvalidOptions("TFIM fidelity scan was empty".into()))?;
    Ok(observation()
        .metric("residual", residual)
        .metric("minimum_fidelity", *minimum_fidelity)
        .metric(
            "minimum_fidelity_interval",
            usize_as_f64(minimum_index + 1)?,
        ))
}

fn pxp_revival() -> Result<Observation, QuSpinError> {
    let sites = 24;
    let basis = UserBasis::builder(sites)
        .states(periodic_blockade_states(sites))
        .operator('x', |state, site| {
            Ok(Some((state ^ (1_u128 << site), Complex64::new(1.0, 0.0))))
        })
        .build()?;
    require_dimension(basis.len(), 103_682)?;
    let hamiltonian = OperatorBuilder::on(&basis)
        .term(OperatorTerm::new(
            "x",
            (0..sites).map(|site| Coupling::new(1.0, vec![site])),
        )?)
        .build(PublicMatrixFormat::Csc)?;
    let neel = (0..sites)
        .step_by(2)
        .fold(0_u128, |state, site| state | (1_u128 << site));
    let neel_index = basis.index(neel)?;
    let mut initial = vec![c(0.0); basis.len()];
    initial[neel_index] = c(1.0);
    let trajectory = evolve(
        &hamiltonian,
        &initial,
        EvolutionOptions {
            times: vec![0.0, 2.4, 4.8, 7.2, 9.6],
            krylov_dimension: 100,
            tolerance: 1.0e-9,
            max_substeps: 100,
            hamiltonian: true,
        },
    )?;
    let norm_error = trajectory
        .states
        .iter()
        .map(|state| (state_norm(state) - 1.0).abs())
        .fold(0.0_f64, f64::max);
    let fidelities: Vec<_> = trajectory
        .states
        .iter()
        .map(|state| state[neel_index].norm_sqr())
        .collect();
    Ok(observation()
        .metric("basis_dimension", usize_as_f64(basis.len())?)
        .metric("norm_error", norm_error)
        .metric("revival_gain", fidelities[2] - fidelities[1]))
}

fn bose_hubbard_mott_quench() -> Result<Observation, QuSpinError> {
    let sites = 11;
    let basis = BosonBasis1D::builder(sites, 3).particles(sites).build()?;
    require_dimension(basis.len(), 25_653)?;
    let bonds = 0..(sites - 1);
    let hamiltonian = OperatorBuilder::on(&basis)
        .terms([
            OperatorTerm::new(
                "+-",
                bonds
                    .clone()
                    .map(|site| Coupling::new(-0.1, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "-+",
                bonds
                    .clone()
                    .map(|site| Coupling::new(-0.1, vec![site, site + 1])),
            )?,
            OperatorTerm::new(
                "nn",
                (0..sites).map(|site| Coupling::new(0.5, vec![site, site])),
            )?,
            OperatorTerm::new("n", (0..sites).map(|site| Coupling::new(-0.5, vec![site])))?,
        ])
        .build(PublicMatrixFormat::Csc)?;
    let mut mott = 0_u128;
    let mut place = 1_u128;
    for _ in 0..sites {
        mott += place;
        place *= 3;
    }
    let mott_index = basis.index(mott)?;
    let mut initial = vec![c(0.0); basis.len()];
    initial[mott_index] = c(1.0);
    let trajectory = evolve(
        &hamiltonian,
        &initial,
        EvolutionOptions {
            times: vec![0.0, 25.0, 50.0, 100.0, 200.0],
            krylov_dimension: 100,
            tolerance: 1.0e-9,
            max_substeps: 1_000,
            hamiltonian: true,
        },
    )?;
    let norm_error = trajectory
        .states
        .iter()
        .map(|state| (state_norm(state) - 1.0).abs())
        .fold(0.0_f64, f64::max);
    let minimum_return = trajectory.states[1..]
        .iter()
        .map(|state| state[mott_index].norm_sqr())
        .fold(1.0_f64, f64::min);
    Ok(observation()
        .metric("norm_error", norm_error)
        .metric("minimum_return_after_t0", minimum_return))
}

fn spinful_kinetic_terms(sites: usize) -> Result<[OperatorTerm; 4], QuSpinError> {
    let bonds = 0..(sites - 1);
    Ok([
        OperatorTerm::new(
            "+-|",
            bonds
                .clone()
                .map(|site| Coupling::new(-1.0, vec![site, site + 1])),
        )?,
        OperatorTerm::new(
            "-+|",
            bonds
                .clone()
                .map(|site| Coupling::new(1.0, vec![site, site + 1])),
        )?,
        OperatorTerm::new(
            "|+-",
            bonds
                .clone()
                .map(|site| Coupling::new(-1.0, vec![site, site + 1])),
        )?,
        OperatorTerm::new(
            "|-+",
            bonds.map(|site| Coupling::new(1.0, vec![site, site + 1])),
        )?,
    ])
}

fn hubbard_interaction(sites: usize, strength: f64) -> Result<OperatorTerm, QuSpinError> {
    OperatorTerm::new(
        "n|n",
        (0..sites).map(|site| Coupling::new(strength, vec![site, site])),
    )
}

fn spinful_hubbard_current_quench() -> Result<Observation, QuSpinError> {
    let sites = 10;
    let basis = SpinfulFermionBasis1D::builder(sites)
        .particles(5, 5)
        .build()?;
    require_dimension(basis.len(), 63_504)?;
    let mut biased_terms = spinful_kinetic_terms(sites)?.to_vec();
    biased_terms.push(hubbard_interaction(sites, 8.0)?);
    biased_terms.extend([
        OperatorTerm::new(
            "n|",
            (0..sites)
                .map(|site| Coupling::new(if site < sites / 2 { -1.5 } else { 1.5 }, vec![site])),
        )?,
        OperatorTerm::new(
            "|n",
            (0..sites)
                .map(|site| Coupling::new(if site < sites / 2 { -1.5 } else { 1.5 }, vec![site])),
        )?,
    ]);
    let biased = OperatorBuilder::on(&basis)
        .terms(biased_terms)
        .build(PublicMatrixFormat::Csc)?;
    let ground = eigsh(
        &biased,
        EigshOptions {
            eigenpairs: 1,
            target: SpectrumTarget::SmallestAlgebraic,
            krylov_dimension: Some(200),
            tolerance: 1.0e-9,
            max_iterations: 240,
            seed: 53,
        },
    )?;
    let mut unbiased_terms = spinful_kinetic_terms(sites)?.to_vec();
    unbiased_terms.push(hubbard_interaction(sites, 8.0)?);
    let unbiased = OperatorBuilder::on(&basis)
        .terms(unbiased_terms)
        .build(PublicMatrixFormat::Csc)?;
    let center = sites / 2 - 1;
    let minus_i = Complex64::new(0.0, -1.0);
    let current = OperatorBuilder::on(&basis)
        .terms([
            OperatorTerm::new("+-|", [Coupling::new(minus_i, vec![center, center + 1])])?,
            OperatorTerm::new("-+|", [Coupling::new(minus_i, vec![center, center + 1])])?,
            OperatorTerm::new("|+-", [Coupling::new(minus_i, vec![center, center + 1])])?,
            OperatorTerm::new("|-+", [Coupling::new(minus_i, vec![center, center + 1])])?,
        ])
        .build(PublicMatrixFormat::Csc)?;
    let trajectory = evolve(
        &unbiased,
        &ground.eigenvectors[0],
        EvolutionOptions {
            times: vec![0.0, 0.5, 1.0, 1.5, 2.0],
            krylov_dimension: 100,
            tolerance: 1.0e-9,
            max_substeps: 100,
            hamiltonian: true,
        },
    )?;
    let mut maximum_current = 0.0_f64;
    let mut applied = vec![c(0.0); basis.len()];
    for state in &trajectory.states {
        current.apply(state, &mut applied)?;
        let expectation: Complex64 = state
            .iter()
            .zip(&applied)
            .map(|(left, right)| left.conj() * *right)
            .sum();
        maximum_current = maximum_current.max(expectation.re.abs());
    }
    Ok(observation()
        .metric("residual", ground.residuals[0])
        .metric("maximum_absolute_current", maximum_current))
}

fn conb_dynamical_structure_factor() -> Result<Observation, QuSpinError> {
    let sites = 16;
    let basis = SpinBasis1D::builder(sites).pauli(false).build()?;
    require_dimension(basis.len(), 65_536)?;
    let transverse_field = 3.21 * 0.057_883_8 * 7.0 / 2.88;
    let bonds = |distance: usize, coefficient: f64, operator: &str| {
        OperatorTerm::new(
            operator,
            (0..sites)
                .map(move |site| Coupling::new(coefficient, vec![site, (site + distance) % sites])),
        )
    };
    let hamiltonian = OperatorBuilder::on(&basis)
        .terms([
            bonds(1, -1.0, "zz")?,
            bonds(1, -0.205, "xx")?,
            bonds(1, -0.205, "yy")?,
            bonds(2, 0.135, "zz")?,
            bonds(2, 0.003, "xx")?,
            bonds(2, 0.003, "yy")?,
            OperatorTerm::new(
                "x",
                (0..sites).map(|site| Coupling::new(-transverse_field, vec![site])),
            )?,
        ])
        .build(PublicMatrixFormat::Csc)?;
    let ground = eigsh(
        &hamiltonian,
        EigshOptions {
            eigenpairs: 1,
            target: SpectrumTarget::SmallestAlgebraic,
            krylov_dimension: Some(180),
            tolerance: 1.0e-9,
            max_iterations: 220,
            seed: 59,
        },
    )?;
    let spin_q = OperatorBuilder::on(&basis)
        .term(OperatorTerm::new(
            "z",
            (0..sites)
                .map(|site| Coupling::new(if site % 2 == 0 { 1.0 } else { -1.0 }, vec![site])),
        )?)
        .build(PublicMatrixFormat::Csc)?;
    let spectrum = spectral_function(
        &hamiltonian,
        &ground.eigenvectors[0],
        &spin_q,
        SpectrumOptions {
            frequencies: (0..=80)
                .map(|index| 4.0 * f64::from(index) / 80.0)
                .collect(),
            reference_energy: ground.eigenvalues[0],
            broadening: 0.05,
            krylov_dimension: 100,
            tolerance: 1.0e-9,
        },
    )?;
    let finite_fraction = usize_as_f64(spectrum.iter().filter(|value| value.is_finite()).count())?
        / usize_as_f64(spectrum.len())?;
    let minimum = spectrum.iter().copied().fold(f64::INFINITY, f64::min);
    let maximum = spectrum.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Ok(observation()
        .metric("residual", ground.residuals[0])
        .metric("minimum_spectrum", minimum)
        .metric("maximum_spectrum", maximum)
        .metric("finite_fraction", finite_fraction))
}

fn triangular_bonds(width: usize, height: usize) -> BTreeSet<(usize, usize)> {
    let site = |x: usize, y: usize| (y % height) * width + (x % width);
    let mut bonds = BTreeSet::new();
    for y in 0..height {
        for x in 0..width {
            let origin = site(x, y);
            for neighbor in [site(x + 1, y), site(x, y + 1), site(x + 1, y + 1)] {
                bonds.insert((origin.min(neighbor), origin.max(neighbor)));
            }
        }
    }
    bonds
}

fn particle_addition_spectrum() -> Result<Observation, QuSpinError> {
    let (width, height) = (6, 3);
    let sites = width * height;
    let source_basis = SpinlessFermionBasis1D::builder(sites)
        .particles(6)
        .build()?;
    let target_basis = SpinlessFermionBasis1D::builder(sites)
        .particles(7)
        .build()?;
    require_dimension(source_basis.len(), 18_564)?;
    require_dimension(target_basis.len(), 31_824)?;
    let bonds = triangular_bonds(width, height);
    if bonds.len() != 54 {
        return Err(QuSpinError::InvalidSector(format!(
            "triangular workflow expected 54 bonds, got {}",
            bonds.len()
        )));
    }
    let hamiltonian_terms = || {
        Ok::<_, QuSpinError>([
            OperatorTerm::new(
                "+-",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(-1.0, vec![left, right])),
            )?,
            OperatorTerm::new(
                "-+",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(1.0, vec![left, right])),
            )?,
            OperatorTerm::new(
                "nn",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(2.0, vec![left, right])),
            )?,
        ])
    };
    let source_hamiltonian = OperatorBuilder::on(&source_basis)
        .terms(hamiltonian_terms()?)
        .build(PublicMatrixFormat::Csc)?;
    let target_hamiltonian = OperatorBuilder::on(&target_basis)
        .terms(hamiltonian_terms()?)
        .build(PublicMatrixFormat::Csc)?;
    let ground = eigsh(
        &source_hamiltonian,
        EigshOptions {
            eigenpairs: 1,
            target: SpectrumTarget::SmallestAlgebraic,
            krylov_dimension: Some(200),
            tolerance: 1.0e-9,
            max_iterations: 240,
            seed: 61,
        },
    )?;
    let probe = OperatorBuilder::between(&source_basis, &target_basis)
        .term(OperatorTerm::new(
            "+",
            [Coupling::new(1.0, vec![sites / 2])],
        )?)
        .build(PublicMatrixFormat::Csc)?;
    let mut transition = vec![c(0.0); target_basis.len()];
    probe.apply(&ground.eigenvectors[0], &mut transition)?;
    let transition_norm = state_norm(&transition);
    let spectrum = spectral_function(
        &target_hamiltonian,
        &ground.eigenvectors[0],
        &probe,
        SpectrumOptions {
            frequencies: (0..=80)
                .map(|index| -4.0 + 16.0 * f64::from(index) / 80.0)
                .collect(),
            reference_energy: ground.eigenvalues[0],
            broadening: 0.1,
            krylov_dimension: 100,
            tolerance: 1.0e-9,
        },
    )?;
    let finite_fraction = usize_as_f64(spectrum.iter().filter(|value| value.is_finite()).count())?
        / usize_as_f64(spectrum.len())?;
    let maximum = spectrum.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Ok(observation()
        .metric("residual", ground.residuals[0])
        .metric("transition_norm", transition_norm)
        .metric("maximum_spectrum", maximum)
        .metric("finite_fraction", finite_fraction))
}

impl WorkflowBackend for QuSpinAdapter {
    type Error = QuSpinError;

    fn language(&self) -> &'static str {
        "rust"
    }

    fn run(&mut self, case: &WorkflowCase) -> Result<Observation, Self::Error> {
        match case.case_id {
            "paper_mbl_shift_invert_l14" => mbl_shift_invert(),
            "paper_xxz_lanczos_quench_l16" => xxz_lanczos_quench(),
            "paper_floquet_heating_l9" => floquet_heating(),
            "paper_spinful_hubbard_l8" => spinful_hubbard(),
            "paper_interacting_ssh_l16" => interacting_ssh(),
            "paper_translation_xxz_l18" => translation_sector_xxz(),
            "paper_tfim_fidelity_l16" => tfim_fidelity_scan(),
            "paper_pxp_revival_l24" => pxp_revival(),
            "paper_bose_hubbard_quench_l11" => bose_hubbard_mott_quench(),
            "paper_hubbard_current_l10" => spinful_hubbard_current_quench(),
            "paper_conb_dsf_l16" => conb_dynamical_structure_factor(),
            "paper_particle_addition_6x3" => particle_addition_spectrum(),
            identifier => Err(QuSpinError::InvalidOptions(format!(
                "unknown private workflow case {identifier}"
            ))),
        }
    }
}
