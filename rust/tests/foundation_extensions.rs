use std::sync::Arc;

use qmbed::archive::{load_basis_zip, save_basis_zip, BasisArchive};
use qmbed::basis::{ErasedState, U256};
use qmbed::dynamics::{DriveStep, Floquet, FloquetSpectrumOptions};
use qmbed::measure::{entanglement_entropy_sector, partial_trace_sector_state, EntropyOrder};
use qmbed::operator::{MatrixFormat, Operator};
use qmbed::Complex64;

#[test]
fn wide_sector_subsystem_contraction_does_not_enumerate_the_parent_space() {
    let states = [
        U256::zero()
            .with_bit(0, true)
            .unwrap()
            .with_bit(199, true)
            .unwrap(),
        U256::zero()
            .with_bit(1, true)
            .unwrap()
            .with_bit(198, true)
            .unwrap(),
    ];
    let scale = 1.0 / 2.0_f64.sqrt();
    let amplitudes = [Complex64::new(scale, 0.0), Complex64::new(scale, 0.0)];
    let reduced = partial_trace_sector_state(&amplitudes, &states, 200, &[0], &[]).unwrap();
    assert!((reduced[0].re - 0.5).abs() < 1.0e-12);
    assert!((reduced[3].re - 0.5).abs() < 1.0e-12);
    let entropy = entanglement_entropy_sector(
        &amplitudes,
        &states,
        200,
        &[0],
        &[],
        EntropyOrder::VonNeumann,
    )
    .unwrap();
    assert!((entropy - 2.0_f64.ln()).abs() < 1.0e-12);
}

#[test]
fn portable_basis_manifest_round_trips_wide_identifiers() {
    let states = [
        "0",
        "3",
        "1606938044258990275541962092341162602522202993782792835301376",
    ]
    .into_iter()
    .map(|value| ErasedState::from_decimal(256, value).unwrap())
    .collect::<Vec<_>>();
    let mut archive = BasisArchive::new(256, states.clone()).unwrap();
    archive.insert_metadata("kind", "spin").unwrap();
    let path = std::env::temp_dir().join(format!(
        "qmbed-benchmark-basis-{}-{}.npz",
        std::process::id(),
        states.len()
    ));
    save_basis_zip(&path, &archive).unwrap();
    let restored = load_basis_zip(&path).unwrap();
    std::fs::remove_file(path).unwrap();
    assert_eq!(restored.width_bits(), 256);
    assert_eq!(restored.states(), states);
    assert_eq!(restored.metadata_value("kind"), Some("spin"));
}

#[test]
fn selected_floquet_spectrum_stays_matrix_free_above_the_dense_cutoff() {
    let dimension = 129;
    let mut energies = vec![-2.0; dimension];
    energies[..3].copy_from_slice(&[0.31, 0.2, 0.45]);
    let hamiltonian = Operator::from_triplets(
        dimension,
        dimension,
        energies
            .iter()
            .copied()
            .enumerate()
            .map(|(index, energy)| (index, index, Complex64::new(energy, 0.0))),
        MatrixFormat::Csc,
    )
    .unwrap();
    let floquet = Floquet::new([DriveStep::new(Arc::new(hamiltonian), 1.0).unwrap()]).unwrap();
    let selected = floquet
        .selected_eigensystem(
            FloquetSpectrumOptions::new(3, 0.31)
                .with_search_dimension(5)
                .with_krylov_dimension(12)
                .with_tolerance(1.0e-11),
        )
        .unwrap();
    assert!(selected.residuals.iter().all(|residual| *residual < 1.0e-8));
    assert!((selected.quasienergies[0] - 0.31).abs() < 1.0e-8);
}
