use qmbed_benchmark::candidate::QmbedAdapter;
use qmbed_benchmark::{
    AssemblyChecks, BasisHandle, BasisSpec, ComplexCoefficient, Coupling, EigshOptions,
    HamiltonianOptions, MatrixFormat, OperatorTerm, QmbedApi, SpectrumTarget,
};
use quspin::Complex64;

#[test]
fn public_candidate_is_reached_through_the_verification_adapter() {
    let adapter = QmbedAdapter;
    let basis = adapter
        .basis(&BasisSpec::Spin {
            sites: 2,
            spin_twice: 1,
            magnetization: None,
            momentum: None,
            parity: None,
            pauli: false,
        })
        .unwrap();
    assert_eq!(basis.dimension(), 4);
    for index in 0..basis.dimension() {
        let state = basis.state_at(index).unwrap();
        assert_eq!(basis.state_index(state).unwrap(), index);
    }

    let coupling = |coefficient, sites| Coupling {
        coefficient: ComplexCoefficient::from(coefficient),
        sites,
    };
    let terms = [
        OperatorTerm {
            operator: "zz".into(),
            couplings: vec![coupling(1.0, vec![0, 1])],
        },
        OperatorTerm {
            operator: "+-".into(),
            couplings: vec![coupling(0.5, vec![0, 1])],
        },
        OperatorTerm {
            operator: "-+".into(),
            couplings: vec![coupling(0.5, vec![0, 1])],
        },
    ];
    let hamiltonian = adapter
        .hamiltonian(
            &basis,
            &terms,
            HamiltonianOptions {
                format: MatrixFormat::Csc,
                checks: AssemblyChecks::default(),
            },
        )
        .unwrap();
    let result = adapter
        .eigsh(
            &hamiltonian,
            None::<&Vec<Complex64>>,
            &EigshOptions {
                eigenpairs: 2,
                target: SpectrumTarget::SmallestAlgebraic,
                krylov_dimension: None,
                tolerance: 1.0e-12,
                max_iterations: 100,
            },
        )
        .unwrap();
    assert!((result.eigenvalues[0] + 0.75).abs() < 1.0e-12);
    assert!((result.eigenvalues[1] - 0.25).abs() < 1.0e-12);
}
