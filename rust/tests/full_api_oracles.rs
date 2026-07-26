use std::collections::HashMap;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use qmbed::archive::{load_zip, save_zip, OperatorArchive};
use qmbed::basis::{
    Basis, BasisProjector, BosonBasis1D, ClosureSymmetryMap, GeneralBasis, PhotonBasis,
    SpinBasis1D, SpinfulFermionBasis1D, SymmetryReducer, UserBasis, WideSpinBasis256, U256,
};
use qmbed::measure::{
    diagonal_ensemble_density, entanglement_entropy, entanglement_spectrum_density,
    partial_trace_density, EntropyOrder,
};
use qmbed::operator::{
    AssemblyChecks, Coupling, LinearOperator, MatrixFormat, Operator, OperatorBuilder,
    OperatorSpec, QuantumComponent, QuantumLinearOperator, QuantumOperator, TimeOperator,
};
use qmbed::solve::{
    eigsh, lanczos_full, EigshOptions, ExpmMultiplyParallel, LanczosOptions, SpectrumTarget,
};
use qmbed::Complex64;

fn close(actual: Complex64, expected: Complex64, tolerance: f64) {
    assert!(
        (actual - expected).norm() <= tolerance,
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

fn inner(left: &[Complex64], right: &[Complex64]) -> Complex64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| left.conj() * *right)
        .sum()
}

fn small_u128(value: u128) -> f64 {
    f64::from(u32::try_from(value).unwrap())
}

fn small_usize(value: usize) -> f64 {
    f64::from(u32::try_from(value).unwrap())
}

fn unchecked() -> AssemblyChecks {
    AssemblyChecks {
        hermiticity: false,
        particle_conservation: false,
        symmetry_compatibility: true,
    }
}

#[test]
fn held_out_general_symmetry_projectors_are_isometric_and_orthogonal() {
    let reflection = || {
        ClosureSymmetryMap::new(2, |state: u128| {
            let mut reflected = 0_u128;
            for site in 0..5 {
                reflected |= ((state >> site) & 1) << (4 - site);
            }
            Ok((reflected, Complex64::new(1.0, 0.0)))
        })
        .unwrap()
    };
    let even = GeneralBasis::new(
        SpinBasis1D::builder(5).build().unwrap(),
        SymmetryReducer::new().with_map(reflection(), 0),
    )
    .unwrap();
    let odd = GeneralBasis::new(
        SpinBasis1D::builder(5).build().unwrap(),
        SymmetryReducer::new().with_map(reflection(), 1),
    )
    .unwrap();
    assert_eq!(even.len() + odd.len(), 32);
    let even_projector = BasisProjector::from_general(&even).unwrap();
    let odd_projector = BasisProjector::from_general(&odd).unwrap();
    for column in 0..even.len() {
        let mut vector = vec![Complex64::new(0.0, 0.0); even.len()];
        vector[column] = Complex64::new(1.0, 0.0);
        let lifted = even_projector.lifted(&vector).unwrap();
        close(inner(&lifted, &lifted), Complex64::new(1.0, 0.0), 2.0e-13);
        for other in 0..odd.len() {
            let mut odd_vector = vec![Complex64::new(0.0, 0.0); odd.len()];
            odd_vector[other] = Complex64::new(1.0, 0.0);
            close(
                inner(&lifted, &odd_projector.lifted(&odd_vector).unwrap()),
                Complex64::new(0.0, 0.0),
                2.0e-13,
            );
        }
    }
}

#[test]
fn held_out_spinful_majorana_and_higher_spin_algebras_close() {
    let union = SpinfulFermionBasis1D::builder(3)
        .particle_sectors([(1, 0), (0, 2), (2, 1)])
        .build()
        .unwrap();
    assert_eq!(union.len(), 3 + 3 + 9);
    let fermions = SpinfulFermionBasis1D::builder(2).build().unwrap();
    let x = OperatorBuilder::on(&fermions)
        .term(OperatorSpec::new("|x", [Coupling::new(1.0, vec![1])]).unwrap())
        .build(MatrixFormat::Csc)
        .unwrap();
    let y = OperatorBuilder::on(&fermions)
        .term(OperatorSpec::new("|y", [Coupling::new(1.0, vec![1])]).unwrap())
        .build(MatrixFormat::Csc)
        .unwrap();
    assert_eq!(
        x.pow(2).unwrap().diagonal(),
        vec![Complex64::new(1.0, 0.0); 16]
    );
    assert_eq!(
        y.pow(2).unwrap().diagonal(),
        vec![Complex64::new(1.0, 0.0); 16]
    );
    assert_eq!(
        x.product(&y)
            .unwrap()
            .add(&y.product(&x).unwrap())
            .unwrap()
            .nnz(),
        0
    );

    let spin = SpinBasis1D::builder(1).spin_twice(3).build().unwrap();
    let plus = OperatorBuilder::on(&spin)
        .term(OperatorSpec::new("+", [Coupling::new(1.0, vec![0])]).unwrap())
        .checks(unchecked())
        .build(MatrixFormat::Dense)
        .unwrap();
    let minus = OperatorBuilder::on(&spin)
        .term(OperatorSpec::new("-", [Coupling::new(1.0, vec![0])]).unwrap())
        .checks(unchecked())
        .build(MatrixFormat::Dense)
        .unwrap();
    let z = OperatorBuilder::on(&spin)
        .term(OperatorSpec::new("z", [Coupling::new(1.0, vec![0])]).unwrap())
        .build(MatrixFormat::Dense)
        .unwrap();
    let commutator = plus
        .product(&minus)
        .unwrap()
        .subtract(&minus.product(&plus).unwrap())
        .unwrap();
    for (actual, expected) in commutator
        .to_dense()
        .iter()
        .zip(z.scaled(2.0).unwrap().to_dense())
    {
        close(*actual, expected, 3.0e-13);
    }
}

#[test]
fn held_out_wide_and_branching_basis_actions_assemble_directly() {
    let wide = WideSpinBasis256::new(220, Some(1), false).unwrap();
    let vacuum = WideSpinBasis256::new(220, Some(0), false).unwrap();
    let lowering = OperatorBuilder::between(&wide, &vacuum)
        .term(OperatorSpec::new("-", [Coupling::new(-0.7, vec![219])]).unwrap())
        .build(MatrixFormat::Csc)
        .unwrap();
    let high = U256::zero().with_bit(219, true).unwrap();
    close(
        lowering.to_dense()[wide.index(high).unwrap()],
        Complex64::new(-0.7, 0.0),
        1.0e-15,
    );

    let basis = UserBasis::builder(1)
        .states([0_u128, 1, 2, 3])
        .branching_operator('m', |state, _| {
            Ok((0..4)
                .filter(|target| (target + state) % 2 == 0)
                .map(|target| {
                    (
                        target,
                        Complex64::new(0.2 * small_u128(target + 1), -0.1 * small_u128(state)),
                    )
                })
                .collect())
        })
        .build()
        .unwrap();
    let operator = OperatorBuilder::on(&basis)
        .term(OperatorSpec::new("m", [Coupling::new(1.0, vec![0])]).unwrap())
        .checks(unchecked())
        .build(MatrixFormat::Csc)
        .unwrap();
    for row in 0..4 {
        for column in 0..4 {
            let expected = if (row + column) % 2 == 0 {
                Complex64::new(0.2 * small_usize(row + 1), -0.1 * small_usize(column))
            } else {
                Complex64::new(0.0, 0.0)
            };
            close(operator.to_dense()[row * 4 + column], expected, 1.0e-15);
        }
    }
}

#[test]
fn held_out_photon_sector_and_exchange_match_the_explicit_two_state_model() {
    let matter = SpinBasis1D::builder(1).build().unwrap();
    let photon = BosonBasis1D::builder(1, 5).build().unwrap();
    let basis = PhotonBasis::fixed_total_excitations(matter, photon, 3, |state| {
        state.count_ones() as usize
    })
    .unwrap();
    assert_eq!(basis.len(), 2);
    let exchange = OperatorBuilder::on(&basis)
        .terms([
            OperatorSpec::new("+|-", [Coupling::new(1.0, vec![0, 0])]).unwrap(),
            OperatorSpec::new("-|+", [Coupling::new(1.0, vec![0, 0])]).unwrap(),
        ])
        .build(MatrixFormat::Dense)
        .unwrap();
    let expected = 3.0_f64.sqrt();
    close(
        exchange.to_dense()[1],
        Complex64::new(expected, 0.0),
        2.0e-14,
    );
    close(
        exchange.to_dense()[2],
        Complex64::new(expected, 0.0),
        2.0e-14,
    );
}

#[test]
fn held_out_dynamic_parameterized_and_linear_operator_semantics_agree() {
    let left = TimeOperator::new((2, 2), |time, input, output| {
        output[0] = time.cos() * input[1];
        output[1] = time.cos() * input[0];
        Ok(())
    })
    .unwrap();
    let right = TimeOperator::new((2, 2), |time, input, output| {
        output[0] = time * input[0];
        output[1] = -time * input[1];
        Ok(())
    })
    .unwrap();
    let time = 0.37;
    let commutator = left.commutator(&right).unwrap();
    let dense_left = left.evaluate(time, MatrixFormat::Dense).unwrap();
    let dense_right = right.evaluate(time, MatrixFormat::Dense).unwrap();
    let expected = dense_left
        .product(&dense_right)
        .unwrap()
        .subtract(&dense_right.product(&dense_left).unwrap())
        .unwrap();
    for (actual, expected) in commutator
        .evaluate(time, MatrixFormat::Dense)
        .unwrap()
        .to_dense()
        .iter()
        .zip(expected.to_dense())
    {
        close(*actual, expected, 2.0e-15);
    }

    let x = dense_left;
    let z = dense_right;
    let parameterized = QuantumOperator::new([
        QuantumComponent::parameter("x", x.clone()),
        QuantumComponent::parameter("z", z.clone()),
    ])
    .unwrap();
    let evaluated = parameterized
        .evaluate(
            &HashMap::from([("x".to_string(), Complex64::new(0.2, 0.0))]),
            MatrixFormat::Dense,
        )
        .unwrap();
    let expected = x.scaled(0.2).unwrap().add(&z).unwrap();
    assert_eq!(evaluated.to_dense(), expected.to_dense());

    let mut linear =
        QuantumLinearOperator::new(x, vec![Complex64::new(0.3, 0.0), Complex64::new(-0.2, 0.0)])
            .unwrap();
    linear
        .set_diagonal(vec![Complex64::new(-0.4, 0.0), Complex64::new(0.5, 0.0)])
        .unwrap();
    let materialized = linear.materialize(MatrixFormat::Dense).unwrap();
    assert_eq!(
        materialized.diagonal(),
        vec![Complex64::new(-0.4, 0.0), Complex64::new(0.5, 0.0)]
    );
}

#[test]
fn held_out_solver_and_general_exponential_oracles_pass() {
    let diagonal = Operator::from_triplets(
        5,
        5,
        [-4.0, -0.3, 0.2, 2.5, 7.0]
            .into_iter()
            .enumerate()
            .map(|(index, value)| (index, index, Complex64::new(value, 0.0))),
        MatrixFormat::Csc,
    )
    .unwrap();
    let selected = eigsh(
        &diagonal,
        EigshOptions::new(3, SpectrumTarget::SmallestMagnitude)
            .with_tolerance(1.0e-12)
            .with_max_iterations(100)
            .with_seed(19),
    )
    .unwrap();
    for (actual, expected) in selected.eigenvalues.iter().zip([0.2, -0.3, 2.5]) {
        assert!(
            (actual - expected).abs() <= 1.0e-12 * expected.abs().max(1.0),
            "eigenvalue mismatch: actual={actual:.16e}, expected={expected:.16e}"
        );
    }

    let initial = vec![
        Complex64::new(1.0, 0.0),
        Complex64::new(-0.4, 0.2),
        Complex64::new(0.3, -0.7),
        Complex64::new(0.1, 0.5),
        Complex64::new(-0.2, 0.0),
    ];
    let decomposition = lanczos_full(
        &diagonal,
        &initial,
        LanczosOptions::new(5).with_tolerance(1.0e-13),
    )
    .unwrap();
    for (index, vector) in decomposition.basis.iter().enumerate() {
        for (other, other_vector) in decomposition.basis.iter().enumerate() {
            close(
                inner(vector, other_vector),
                Complex64::new(if index == other { 1.0 } else { 0.0 }, 0.0),
                3.0e-12,
            );
        }
    }

    let rotation = Arc::new(
        Operator::from_dense(
            2,
            2,
            vec![
                Complex64::new(0.0, 0.0),
                Complex64::new(-1.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
            ],
        )
        .unwrap(),
    );
    let plan =
        ExpmMultiplyParallel::new(rotation, Complex64::new(0.43, 0.0), 32, 1.0e-14, 100).unwrap();
    let mut output = vec![Complex64::new(0.0, 0.0); 2];
    plan.apply(
        &[Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        &mut output,
    )
    .unwrap();
    close(output[0], Complex64::new(0.43_f64.cos(), 0.0), 2.0e-13);
    close(output[1], Complex64::new(0.43_f64.sin(), 0.0), 2.0e-13);
}

#[test]
fn held_out_mixed_measurements_and_named_archives_round_trip() {
    let amplitude = 1.0 / 2.0_f64.sqrt();
    let bell = [
        Complex64::new(amplitude, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(amplitude, 0.0),
    ];
    let density: Vec<_> = bell
        .iter()
        .flat_map(|left| bell.iter().map(move |right| *left * right.conj()))
        .collect();
    let reduced = partial_trace_density(&density, 2, 2).unwrap();
    close(reduced[0], Complex64::new(0.5, 0.0), 2.0e-15);
    close(reduced[3], Complex64::new(0.5, 0.0), 2.0e-15);
    assert_eq!(
        entanglement_spectrum_density(&density, 2, 2).unwrap(),
        vec![0.5, 0.5]
    );
    close(
        Complex64::new(
            entanglement_entropy(&bell, 2, 2, EntropyOrder::VonNeumann).unwrap(),
            0.0,
        ),
        Complex64::new(2.0_f64.ln(), 0.0),
        2.0e-15,
    );
    let eigenvectors = (0..4)
        .map(|column| {
            (0..4)
                .map(|row| Complex64::new(if row == column { 1.0 } else { 0.0 }, 0.0))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let ensemble =
        diagonal_ensemble_density(&[-1.0, -0.2, 0.7, 2.0], &eigenvectors, &density).unwrap();
    close(
        Complex64::new(ensemble.probabilities.iter().sum(), 0.0),
        Complex64::new(1.0, 0.0),
        2.0e-15,
    );

    let dense = Operator::from_dense(2, 2, reduced).unwrap();
    let sparse = Operator::from_triplets(
        2,
        2,
        [
            (0, 1, Complex64::new(0.0, 0.75)),
            (1, 0, Complex64::new(0.0, -0.75)),
        ],
        MatrixFormat::Csr,
    )
    .unwrap();
    let mut archive = OperatorArchive::new();
    archive
        .insert("density", dense.clone(), Some(Complex64::new(1.0, 0.0)))
        .unwrap();
    archive.insert("current", sparse.clone(), None).unwrap();
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("qmbed-benchmark-{nonce}.npz"));
    save_zip(&path, &archive).unwrap();
    let restored = load_zip(&path).unwrap();
    assert_eq!(
        restored.get("density").unwrap().operator.to_dense(),
        dense.to_dense()
    );
    assert_eq!(
        restored.get("current").unwrap().operator.format(),
        MatrixFormat::Csr
    );
    assert_eq!(
        restored.get("current").unwrap().operator.triplets(),
        sparse.triplets()
    );
    std::fs::remove_file(path).unwrap();
}
