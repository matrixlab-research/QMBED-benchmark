use qmbed::basis::{Basis, BosonBasis1D, UserBasis};
use qmbed::operator::{Coupling, MatrixFormat, Operator, OperatorBuilder, OperatorTerm};
use qmbed::solve::{evolve_with_diagnostics, lanczos_ritz, EvolutionOptions, LanczosOptions};
use qmbed::{Complex64, QmbedError};

struct EvolutionCase {
    name: &'static str,
    operator: Operator,
    initial: Vec<Complex64>,
    times: Vec<f64>,
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

fn pxp_case() -> Result<EvolutionCase, QmbedError> {
    let sites = 24;
    let basis = UserBasis::builder(sites)
        .states(periodic_blockade_states(sites))
        .operator('x', |state, site| {
            Ok(Some((state ^ (1_u128 << site), Complex64::new(1.0, 0.0))))
        })
        .build()?;
    let operator = OperatorBuilder::on(&basis)
        .term(OperatorTerm::new(
            "x",
            (0..sites).map(|site| Coupling::new(1.0, vec![site])),
        )?)
        .build(MatrixFormat::Csc)?;
    let neel = (0..sites)
        .step_by(2)
        .fold(0_u128, |state, site| state | (1_u128 << site));
    let mut initial = vec![Complex64::new(0.0, 0.0); basis.len()];
    initial[basis.index(neel)?] = Complex64::new(1.0, 0.0);
    Ok(EvolutionCase {
        name: "PXP",
        operator,
        initial,
        times: vec![0.0, 2.4, 4.8, 7.2, 9.6],
    })
}

fn bose_hubbard_case() -> Result<EvolutionCase, QmbedError> {
    let sites = 11;
    let basis = BosonBasis1D::builder(sites, 3).particles(sites).build()?;
    let bonds = 0..(sites - 1);
    let operator = OperatorBuilder::on(&basis)
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
        .build(MatrixFormat::Csc)?;
    let mut mott = 0_u128;
    let mut place = 1_u128;
    for _ in 0..sites {
        mott += place;
        place *= 3;
    }
    let mut initial = vec![Complex64::new(0.0, 0.0); basis.len()];
    initial[basis.index(mott)?] = Complex64::new(1.0, 0.0);
    Ok(EvolutionCase {
        name: "Bose-Hubbard",
        operator,
        initial,
        times: vec![0.0, 25.0, 50.0, 100.0, 200.0],
    })
}

fn maximum_state_error(left: &[Vec<Complex64>], right: &[Vec<Complex64>]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left_state, right_state)| {
            left_state
                .iter()
                .zip(right_state)
                .map(|(left_value, right_value)| (*left_value - *right_value).norm_sqr())
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0_f64, f64::max)
}

#[test]
#[ignore = "paper-scale accuracy oracle; exercised in release mode"]
fn adaptive_krylov_repairs_the_legacy_single_projection_error() {
    for case in [pxp_case().unwrap(), bose_hubbard_case().unwrap()] {
        let (candidate, diagnostics) = evolve_with_diagnostics(
            &case.operator,
            &case.initial,
            EvolutionOptions {
                times: case.times.clone(),
                krylov_dimension: 100,
                tolerance: 1.0e-9,
                max_substeps: 10_000,
                hamiltonian: true,
            },
        )
        .unwrap();
        let (reference, _) = evolve_with_diagnostics(
            &case.operator,
            &case.initial,
            EvolutionOptions {
                times: case.times.clone(),
                krylov_dimension: 100,
                tolerance: 1.0e-12,
                max_substeps: 10_000,
                hamiltonian: true,
            },
        )
        .unwrap();
        let legacy_projection = lanczos_ritz(
            &case.operator,
            &case.initial,
            LanczosOptions {
                krylov_dimension: 100,
                tolerance: 1.0e-14,
            },
        )
        .unwrap();
        let legacy = case
            .times
            .iter()
            .map(|time| {
                legacy_projection
                    .exponential_action(Complex64::new(0.0, -*time))
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let candidate_error = maximum_state_error(&candidate.states, &reference.states);
        let legacy_error = maximum_state_error(&legacy, &reference.states);
        println!(
            "{},candidate_error={candidate_error:.12e},legacy_error={legacy_error:.12e},\
             projections={},matvecs={},real_projections={},real_matvecs={}",
            case.name,
            diagnostics.lanczos_projections,
            diagnostics.matrix_vector_products,
            diagnostics.real_lanczos_projections,
            diagnostics.real_matrix_vector_products
        );
        assert!(candidate_error <= 2.0e-8);
        assert!(legacy_error >= 0.1);
        assert!(diagnostics.lanczos_projections > 1);
        assert_eq!(diagnostics.real_lanczos_projections, 1);
        assert_eq!(
            diagnostics.matrix_vector_products,
            diagnostics.lanczos_projections * 100
        );
        assert_eq!(diagnostics.real_matrix_vector_products, 100);
    }
}
