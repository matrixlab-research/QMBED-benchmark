//! Parameterized, paper-shaped workflows used to validate native AD.
//!
//! Each case builds a genuine Hilbert space and Hermitian operator family.
//! The same cases drive analytic Hellmann--Feynman gradients and central
//! finite-difference oracles; no case-specific derivative is implemented here.

use std::collections::BTreeSet;
use std::time::{Duration, Instant};

use qmbed::ad::{
    ground_state_energy_gradient, GradientStatus, GroundStateEnergyGradient, ParameterValues,
};
use qmbed::basis::{
    Basis, BosonBasis1D, SpinBasis1D, SpinfulFermionBasis1D, SpinlessFermionBasis1D, UserBasis,
};
use qmbed::operator::{
    Coupling, MatrixFormat, OperatorBuilder, OperatorSpec, QuantumComponent, QuantumOperator,
};
use qmbed::solve::{eigsh_with_workspace, EigshOptions, EigshWorkspace};
use qmbed::{Complex64, QmbedError};

/// Size profile for pull-request gates and paper-scale timing runs.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AdScale {
    /// Fast enough for every pull request while retaining nontrivial sectors.
    Ci,
    /// Medium ED sizes representative of the repository's paper workflows.
    Paper,
}

/// One parameterized end-to-end AD acceptance workflow.
#[derive(Clone, Debug)]
pub struct AdWorkflow {
    /// Stable CSV and test identifier.
    pub case_id: &'static str,
    /// Scientific model class.
    pub family_name: &'static str,
    /// Paper workflow from which the differentiable task was derived.
    pub source_workflow: &'static str,
    /// Hermitian operator components in stable parameter order.
    pub operator: QuantumOperator,
    /// Evaluation point.
    pub parameters: ParameterValues,
    /// Selected-spectrum controls shared by AD and finite differences.
    pub options: EigshOptions,
    /// Central-difference displacement.
    pub finite_difference_step: f64,
}

/// Correctness and performance evidence for one differentiable workflow.
#[derive(Clone, Debug)]
pub struct AdBenchmarkRow {
    /// Stable workflow identifier.
    pub case_id: &'static str,
    /// Scientific model class.
    pub family_name: &'static str,
    /// Hilbert-space dimension.
    pub dimension: usize,
    /// Number of differentiated parameters.
    pub parameters: usize,
    /// Maximum absolute difference from central finite differences.
    pub maximum_absolute_error: f64,
    /// Maximum componentwise relative difference.
    pub maximum_relative_error: f64,
    /// Lowest spectral separation returned by the AD rule.
    pub spectral_gap: f64,
    /// Maximum ground-state residual.
    pub residual: f64,
    /// Median analytic-gradient wall time.
    pub analytic_seconds: f64,
    /// Median central-finite-difference wall time.
    pub finite_difference_seconds: f64,
    /// `finite_difference_seconds / analytic_seconds`.
    pub speedup: f64,
    /// Number of eigensolves in the analytic rule.
    pub analytic_eigensolves: usize,
    /// Number of eigensolves in the finite-difference oracle.
    pub finite_difference_eigensolves: usize,
}

impl AdBenchmarkRow {
    /// Stable CSV header used by CI artifacts and plot scripts.
    pub const CSV_HEADER: &'static str = "case_id,family,dimension,parameters,max_abs_error,max_rel_error,spectral_gap,residual,analytic_seconds,finite_difference_seconds,speedup,analytic_eigensolves,finite_difference_eigensolves";

    /// Serialize without locale-dependent formatting.
    #[must_use]
    pub fn to_csv_record(&self) -> String {
        format!(
            "{},{},{},{},{:.12e},{:.12e},{:.12e},{:.12e},{:.9},{:.9},{:.6},{},{}",
            self.case_id,
            self.family_name,
            self.dimension,
            self.parameters,
            self.maximum_absolute_error,
            self.maximum_relative_error,
            self.spectral_gap,
            self.residual,
            self.analytic_seconds,
            self.finite_difference_seconds,
            self.speedup,
            self.analytic_eigensolves,
            self.finite_difference_eigensolves
        )
    }
}

fn analytic_gradient(case: &AdWorkflow) -> Result<GroundStateEnergyGradient, QmbedError> {
    ground_state_energy_gradient(
        &case.operator,
        &case.parameters,
        case.options.clone(),
        &mut EigshWorkspace::new(),
    )
}

fn finite_difference_gradient(case: &AdWorkflow) -> Result<Vec<f64>, QmbedError> {
    let mut plan = case.operator.plan(MatrixFormat::Csc)?;
    let mut workspace = EigshWorkspace::new();
    let mut gradient = Vec::with_capacity(case.parameters.values().len());
    for parameter in 0..case.parameters.values().len() {
        let mut positive = case.parameters.values().to_vec();
        let mut negative = case.parameters.values().to_vec();
        positive[parameter].re += case.finite_difference_step;
        negative[parameter].re -= case.finite_difference_step;
        let positive_energy = eigsh_with_workspace(
            plan.evaluate_coefficients(&positive)?,
            case.options.clone(),
            &mut workspace,
        )?
        .eigenvalues[0];
        let negative_energy = eigsh_with_workspace(
            plan.evaluate_coefficients(&negative)?,
            case.options.clone(),
            &mut workspace,
        )?
        .eigenvalues[0];
        gradient.push((positive_energy - negative_energy) / (2.0 * case.finite_difference_step));
    }
    Ok(gradient)
}

fn median(samples: &mut [Duration]) -> f64 {
    samples.sort_unstable();
    samples[samples.len() / 2].as_secs_f64()
}

/// Validate and time one workflow against a central-difference oracle.
///
/// # Errors
///
/// Returns an error when solver execution, derivative validation, or the
/// requested timing protocol fails.
pub fn benchmark_ad_workflow(
    case: &AdWorkflow,
    warmups: usize,
    samples: usize,
) -> Result<AdBenchmarkRow, QmbedError> {
    if samples == 0 {
        return Err(QmbedError::InvalidOptions(
            "AD benchmark requires at least one timing sample".into(),
        ));
    }
    let analytic = analytic_gradient(case)?;
    if analytic.diagnostics.status != GradientStatus::Reliable {
        return Err(QmbedError::InvalidOptions(format!(
            "{} gradient is not reliable: {:?}",
            case.case_id, analytic.diagnostics
        )));
    }
    let finite_difference = finite_difference_gradient(case)?;
    let (maximum_absolute_error, maximum_relative_error) = analytic
        .gradient
        .values()
        .iter()
        .zip(&finite_difference)
        .fold((0.0_f64, 0.0_f64), |(absolute, relative), (left, right)| {
            let difference = (left.re - right).abs();
            (
                absolute.max(difference),
                relative.max(difference / left.re.abs().max(right.abs()).max(1.0e-12)),
            )
        });
    if maximum_absolute_error > 2.0e-5 && maximum_relative_error > 2.0e-4 {
        return Err(QmbedError::InvalidOptions(format!(
            "{} AD/finite-difference mismatch: abs={maximum_absolute_error:e}, rel={maximum_relative_error:e}",
            case.case_id
        )));
    }

    for _ in 0..warmups {
        let _ = analytic_gradient(case)?;
        let _ = finite_difference_gradient(case)?;
    }
    let mut analytic_times = Vec::with_capacity(samples);
    let mut finite_difference_times = Vec::with_capacity(samples);
    for sample in 0..samples {
        if sample.is_multiple_of(2) {
            let start = Instant::now();
            let _ = analytic_gradient(case)?;
            analytic_times.push(start.elapsed());
            let start = Instant::now();
            let _ = finite_difference_gradient(case)?;
            finite_difference_times.push(start.elapsed());
        } else {
            let start = Instant::now();
            let _ = finite_difference_gradient(case)?;
            finite_difference_times.push(start.elapsed());
            let start = Instant::now();
            let _ = analytic_gradient(case)?;
            analytic_times.push(start.elapsed());
        }
    }
    let analytic_seconds = median(&mut analytic_times);
    let finite_difference_seconds = median(&mut finite_difference_times);
    Ok(AdBenchmarkRow {
        case_id: case.case_id,
        family_name: case.family_name,
        dimension: case.operator.shape().0,
        parameters: case.parameters.values().len(),
        maximum_absolute_error,
        maximum_relative_error,
        spectral_gap: analytic.diagnostics.spectral_gap.unwrap_or(f64::NAN),
        residual: analytic.diagnostics.primal_residual.unwrap_or(f64::NAN),
        analytic_seconds,
        finite_difference_seconds,
        speedup: finite_difference_seconds / analytic_seconds,
        analytic_eigensolves: 1,
        finite_difference_eigensolves: 2 * case.parameters.values().len(),
    })
}

fn component<B>(
    basis: &B,
    name: &'static str,
    terms: impl IntoIterator<Item = OperatorSpec>,
) -> Result<QuantumComponent, QmbedError>
where
    B: Basis,
{
    Ok(QuantumComponent::required(
        name,
        OperatorBuilder::on(basis)
            .terms(terms)
            .build(MatrixFormat::Csc)?,
    ))
}

fn workflow(
    case_id: &'static str,
    family_name: &'static str,
    source_workflow: &'static str,
    components: impl IntoIterator<Item = QuantumComponent>,
    parameters: impl IntoIterator<Item = f64>,
) -> Result<AdWorkflow, QmbedError> {
    let operator = QuantumOperator::new(components)?;
    let parameters = ParameterValues::real(&operator, parameters)?;
    let options = EigshOptions::smallest_algebraic(2)
        .with_tolerance(1.0e-9)
        .with_max_iterations(2_000)
        .with_seed(73);
    Ok(AdWorkflow {
        case_id,
        family_name,
        source_workflow,
        operator,
        parameters,
        options,
        finite_difference_step: 1.0e-5,
    })
}

fn nearest_bonds(sites: usize) -> impl Iterator<Item = (usize, usize)> + Clone {
    (0..sites - 1).map(|site| (site, site + 1))
}

fn periodic_bonds(sites: usize) -> impl Iterator<Item = (usize, usize)> + Clone {
    (0..sites).map(move |site| (site, (site + 1) % sites))
}

fn usize_to_f64(value: usize) -> Result<f64, QmbedError> {
    u32::try_from(value).map(f64::from).map_err(|_| {
        QmbedError::InvalidOptions(format!(
            "benchmark integer {value} exceeds exact f64 conversion range"
        ))
    })
}

fn spin_tfim(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 8 } else { 14 };
    let basis = SpinBasis1D::builder(sites).pauli(true).build()?;
    workflow(
        "ad_tfim_ground_energy",
        "spin / transverse-field Ising",
        "paper_tfim_fidelity_l16",
        [
            component(
                &basis,
                "exchange",
                [OperatorSpec::new(
                    "zz",
                    nearest_bonds(sites)
                        .map(|(left, right)| Coupling::new(-1.0, vec![left, right])),
                )?],
            )?,
            component(
                &basis,
                "field",
                [OperatorSpec::new(
                    "x",
                    (0..sites).map(|site| Coupling::new(-1.0, vec![site])),
                )?],
            )?,
            component(
                &basis,
                "longitudinal",
                [OperatorSpec::new(
                    "z",
                    (0..sites).map(|site| Coupling::new(-1.0, vec![site])),
                )?],
            )?,
        ],
        [1.0, 0.91, 0.07],
    )
}

fn spin_xxz(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 10 } else { 16 };
    let basis = SpinBasis1D::builder(sites).up(sites / 2).build()?;
    workflow(
        "ad_xxz_sector_ground_energy",
        "spin / fixed-magnetization XXZ",
        "paper_xxz_lanczos_quench_l16",
        [
            component(
                &basis,
                "exchange_xy",
                [
                    OperatorSpec::new(
                        "+-",
                        nearest_bonds(sites)
                            .map(|(left, right)| Coupling::new(0.5, vec![left, right])),
                    )?,
                    OperatorSpec::new(
                        "-+",
                        nearest_bonds(sites)
                            .map(|(left, right)| Coupling::new(0.5, vec![left, right])),
                    )?,
                ],
            )?,
            component(
                &basis,
                "anisotropy",
                [OperatorSpec::new(
                    "zz",
                    nearest_bonds(sites).map(|(left, right)| Coupling::new(1.0, vec![left, right])),
                )?],
            )?,
            component(
                &basis,
                "staggered_field",
                [OperatorSpec::new(
                    "z",
                    (0..sites).map(|site| {
                        Coupling::new(if site.is_multiple_of(2) { 1.0 } else { -1.0 }, vec![site])
                    }),
                )?],
            )?,
        ],
        [1.0, 0.8, 0.13],
    )
}

fn spin_j1_j2(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 10 } else { 16 };
    let basis = SpinBasis1D::builder(sites).up(sites / 2).build()?;
    let heisenberg = |distance: usize| -> Result<Vec<OperatorSpec>, QmbedError> {
        let bonds = (0..sites - distance).map(move |site| (site, site + distance));
        Ok(vec![
            OperatorSpec::new(
                "+-",
                bonds
                    .clone()
                    .map(|(left, right)| Coupling::new(0.5, vec![left, right])),
            )?,
            OperatorSpec::new(
                "-+",
                bonds
                    .clone()
                    .map(|(left, right)| Coupling::new(0.5, vec![left, right])),
            )?,
            OperatorSpec::new(
                "zz",
                bonds.map(|(left, right)| Coupling::new(1.0, vec![left, right])),
            )?,
        ])
    };
    workflow(
        "ad_j1_j2_frustrated_ground_energy",
        "spin / frustrated J1-J2",
        "LKM frustrated-spin workflows",
        [
            component(&basis, "j1", heisenberg(1)?)?,
            component(&basis, "j2", heisenberg(2)?)?,
            component(
                &basis,
                "edge_field",
                [OperatorSpec::new(
                    "z",
                    [
                        Coupling::new(1.0, vec![0]),
                        Coupling::new(-1.0, vec![sites - 1]),
                    ],
                )?],
            )?,
        ],
        [1.0, 0.42, 0.03],
    )
}

fn spin_translation(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 10 } else { 18 };
    let basis = SpinBasis1D::builder(sites)
        .up(sites / 2)
        .momentum(0)
        .build()?;
    workflow(
        "ad_translation_xxz_ground_energy",
        "spin / translation sector",
        "paper_translation_xxz_l18",
        [
            component(
                &basis,
                "exchange_xy",
                [
                    OperatorSpec::new(
                        "+-",
                        periodic_bonds(sites)
                            .map(|(left, right)| Coupling::new(0.5, vec![left, right])),
                    )?,
                    OperatorSpec::new(
                        "-+",
                        periodic_bonds(sites)
                            .map(|(left, right)| Coupling::new(0.5, vec![left, right])),
                    )?,
                ],
            )?,
            component(
                &basis,
                "anisotropy",
                [OperatorSpec::new(
                    "zz",
                    periodic_bonds(sites)
                        .map(|(left, right)| Coupling::new(1.0, vec![left, right])),
                )?],
            )?,
        ],
        [1.0, 0.73],
    )
}

fn spinless_hopping_component<B: Basis<State = u128>>(
    basis: &B,
    name: &'static str,
    bonds: impl IntoIterator<Item = (usize, usize)>,
) -> Result<QuantumComponent, QmbedError> {
    let bonds: Vec<_> = bonds.into_iter().collect();
    component(
        basis,
        name,
        [
            OperatorSpec::new(
                "+-",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(-1.0, vec![left, right])),
            )?,
            OperatorSpec::new(
                "-+",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(1.0, vec![left, right])),
            )?,
        ],
    )
}

fn spinless_ssh(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 10 } else { 16 };
    let basis = SpinlessFermionBasis1D::builder(sites)
        .particles(sites / 2)
        .build()?;
    workflow(
        "ad_interacting_ssh_ground_energy",
        "spinless fermion / interacting SSH",
        "paper_interacting_ssh_l16",
        [
            spinless_hopping_component(
                &basis,
                "even_hopping",
                nearest_bonds(sites).filter(|(left, _)| left.is_multiple_of(2)),
            )?,
            spinless_hopping_component(
                &basis,
                "odd_hopping",
                nearest_bonds(sites).filter(|(left, _)| !left.is_multiple_of(2)),
            )?,
            component(
                &basis,
                "interaction",
                [OperatorSpec::new(
                    "nn",
                    nearest_bonds(sites).map(|(left, right)| Coupling::new(1.0, vec![left, right])),
                )?],
            )?,
        ],
        [0.61, 1.0, 1.7],
    )
}

fn spinless_disorder(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 10 } else { 14 };
    let basis = SpinlessFermionBasis1D::builder(sites)
        .particles(sites / 2)
        .build()?;
    let disorder_couplings = (0..sites)
        .map(|site| {
            let value = usize_to_f64((37 * site + 11) % 29)? / 14.0 - 1.0;
            Ok(Coupling::new(value, vec![site]))
        })
        .collect::<Result<Vec<_>, QmbedError>>()?;
    workflow(
        "ad_disordered_tv_ground_energy",
        "spinless fermion / disordered t-V",
        "paper_mbl_shift_invert_l14",
        [
            spinless_hopping_component(&basis, "hopping", periodic_bonds(sites))?,
            component(
                &basis,
                "interaction",
                [OperatorSpec::new(
                    "nn",
                    periodic_bonds(sites)
                        .map(|(left, right)| Coupling::new(1.0, vec![left, right])),
                )?],
            )?,
            component(
                &basis,
                "disorder",
                [OperatorSpec::new("n", disorder_couplings)?],
            )?,
        ],
        [1.0, 1.4, 2.1],
    )
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

fn spinless_triangular(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let (width, height, particles) = if scale == AdScale::Ci {
        (3, 3, 3)
    } else {
        (4, 3, 4)
    };
    let sites = width * height;
    let bonds = triangular_bonds(width, height);
    let basis = SpinlessFermionBasis1D::builder(sites)
        .particles(particles)
        .build()?;
    workflow(
        "ad_triangular_fermion_ground_energy",
        "spinless fermion / triangular lattice",
        "paper_particle_addition_6x3",
        [
            spinless_hopping_component(&basis, "hopping", bonds.iter().copied())?,
            component(
                &basis,
                "interaction",
                [OperatorSpec::new(
                    "nn",
                    bonds
                        .iter()
                        .map(|&(left, right)| Coupling::new(1.0, vec![left, right])),
                )?],
            )?,
            component(
                &basis,
                "pinning",
                [OperatorSpec::new(
                    "n",
                    [
                        Coupling::new(1.0, vec![0]),
                        Coupling::new(-0.7, vec![sites - 1]),
                    ],
                )?],
            )?,
        ],
        [1.0, 1.8, 0.05],
    )
}

fn spinful_kinetic<B: Basis<State = u128>>(
    basis: &B,
    sites: usize,
) -> Result<QuantumComponent, QmbedError> {
    let bonds: Vec<_> = nearest_bonds(sites).collect();
    component(
        basis,
        "hopping",
        [
            OperatorSpec::new(
                "+-|",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(-1.0, vec![left, right])),
            )?,
            OperatorSpec::new(
                "-+|",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(1.0, vec![left, right])),
            )?,
            OperatorSpec::new(
                "|+-",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(-1.0, vec![left, right])),
            )?,
            OperatorSpec::new(
                "|-+",
                bonds
                    .iter()
                    .map(|&(left, right)| Coupling::new(1.0, vec![left, right])),
            )?,
        ],
    )
}

fn spinful_hubbard(scale: AdScale, ionic: bool) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 4 } else { 8 };
    let basis = SpinfulFermionBasis1D::builder(sites)
        .particles(sites / 2, sites / 2)
        .build()?;
    let identifier = if ionic {
        "ad_ionic_hubbard_ground_energy"
    } else {
        "ad_hubbard_ground_energy"
    };
    let family = if ionic {
        "spinful fermion / ionic Hubbard"
    } else {
        "spinful fermion / Hubbard"
    };
    workflow(
        identifier,
        family,
        "paper_spinful_hubbard_l8",
        [
            spinful_kinetic(&basis, sites)?,
            component(
                &basis,
                "interaction",
                [OperatorSpec::new(
                    "n|n",
                    (0..sites).map(|site| Coupling::new(1.0, vec![site, site])),
                )?],
            )?,
            component(
                &basis,
                "ionic_potential",
                [
                    OperatorSpec::new(
                        "n|",
                        (0..sites).map(|site| {
                            Coupling::new(
                                if site.is_multiple_of(2) { 1.0 } else { -1.0 },
                                vec![site],
                            )
                        }),
                    )?,
                    OperatorSpec::new(
                        "|n",
                        (0..sites).map(|site| {
                            Coupling::new(
                                if site.is_multiple_of(2) { 1.0 } else { -1.0 },
                                vec![site],
                            )
                        }),
                    )?,
                ],
            )?,
        ],
        if ionic {
            [1.0, 4.0, 0.7]
        } else {
            [1.0, 4.0, 0.03]
        },
    )
}

fn bose_hubbard(scale: AdScale, trapped: bool) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 5 } else { 9 };
    let basis = BosonBasis1D::builder(sites, 3).particles(sites).build()?;
    let center = usize_to_f64(sites - 1)? / 2.0;
    let trap_couplings = (0..sites)
        .map(|site| {
            let distance = usize_to_f64(site)? - center;
            Ok(Coupling::new(distance * distance, vec![site]))
        })
        .collect::<Result<Vec<_>, QmbedError>>()?;
    workflow(
        if trapped {
            "ad_trapped_bose_hubbard_ground_energy"
        } else {
            "ad_bose_hubbard_ground_energy"
        },
        if trapped {
            "boson / trapped Bose-Hubbard"
        } else {
            "boson / Bose-Hubbard"
        },
        "paper_bose_hubbard_quench_l11",
        [
            component(
                &basis,
                "hopping",
                [
                    OperatorSpec::new(
                        "+-",
                        nearest_bonds(sites)
                            .map(|(left, right)| Coupling::new(-1.0, vec![left, right])),
                    )?,
                    OperatorSpec::new(
                        "-+",
                        nearest_bonds(sites)
                            .map(|(left, right)| Coupling::new(-1.0, vec![left, right])),
                    )?,
                ],
            )?,
            component(
                &basis,
                "interaction",
                [
                    OperatorSpec::new(
                        "nn",
                        (0..sites).map(|site| Coupling::new(0.5, vec![site, site])),
                    )?,
                    OperatorSpec::new("n", (0..sites).map(|site| Coupling::new(-0.5, vec![site])))?,
                ],
            )?,
            component(&basis, "trap", [OperatorSpec::new("n", trap_couplings)?])?,
        ],
        if trapped {
            [0.16, 1.0, 0.025]
        } else {
            [0.16, 1.0, 0.001]
        },
    )
}

fn periodic_blockade_states(sites: usize) -> Vec<u128> {
    (0..(1_u128 << sites))
        .filter(|state| {
            (0..sites).all(|site| {
                let next = (site + 1) % sites;
                (state & (1_u128 << site)) == 0 || (state & (1_u128 << next)) == 0
            })
        })
        .collect()
}

fn pxp(scale: AdScale) -> Result<AdWorkflow, QmbedError> {
    let sites = if scale == AdScale::Ci { 10 } else { 20 };
    let basis = UserBasis::builder(sites)
        .states(periodic_blockade_states(sites))
        .operator('x', |state, site| {
            Ok(Some((state ^ (1_u128 << site), Complex64::new(1.0, 0.0))))
        })
        .operator('n', |state, site| {
            Ok(Some((
                state,
                Complex64::new(f64::from(u8::from(((state >> site) & 1) != 0)), 0.0),
            )))
        })
        .build()?;
    workflow(
        "ad_pxp_detuning_ground_energy",
        "constrained basis / PXP",
        "paper_pxp_revival_l24",
        [
            component(
                &basis,
                "flip",
                [OperatorSpec::new(
                    "x",
                    (0..sites).map(|site| Coupling::new(1.0, vec![site])),
                )?],
            )?,
            component(
                &basis,
                "detuning",
                [OperatorSpec::new(
                    "n",
                    (0..sites).map(|site| Coupling::new(1.0, vec![site])),
                )?],
            )?,
            component(
                &basis,
                "staggered_detuning",
                [OperatorSpec::new(
                    "n",
                    (0..sites).map(|site| {
                        Coupling::new(if site.is_multiple_of(2) { 1.0 } else { -1.0 }, vec![site])
                    }),
                )?],
            )?,
        ],
        [1.0, 0.31, 0.07],
    )
}

/// Build all twelve AD acceptance workflows.
///
/// # Errors
///
/// Returns an error if a basis, term, symmetry sector, or parameterized
/// Hamiltonian cannot be constructed.
pub fn ad_workflows(scale: AdScale) -> Result<Vec<AdWorkflow>, QmbedError> {
    Ok(vec![
        spin_tfim(scale)?,
        spin_xxz(scale)?,
        spin_j1_j2(scale)?,
        spin_translation(scale)?,
        spinless_ssh(scale)?,
        spinless_disorder(scale)?,
        spinless_triangular(scale)?,
        spinful_hubbard(scale, false)?,
        spinful_hubbard(scale, true)?,
        bose_hubbard(scale, false)?,
        bose_hubbard(scale, true)?,
        pxp(scale)?,
    ])
}
