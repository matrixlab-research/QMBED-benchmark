use num_bigint::BigUint;
use qmbed::basis::{
    state_from_biguint, state_to_biguint, Basis, BasisProjector, BosonBasis1D, ClosureSymmetryMap,
    GeneralBasis, SpinBasis1D, SpinlessFermionBasis1D, SymmetryReducer, WideSpinBasis1024, U1024,
};
use qmbed::measure::{
    entanglement_entropy_density_subsystem, partial_trace_subsystem, EntropyOrder,
};
use qmbed::operator::{
    apply_sector_shift, Coupling, LinearOperator, MatrixFormat, Operator, OperatorBuilder,
    OperatorSpec,
};
use qmbed::Complex64;

fn periodic_heisenberg(sites: usize) -> Vec<OperatorSpec> {
    let bonds = |coefficient| {
        (0..sites)
            .map(|site| Coupling::new(coefficient, vec![site, (site + 1) % sites]))
            .collect::<Vec<_>>()
    };
    vec![
        OperatorSpec::new("zz", bonds(1.0)).unwrap(),
        OperatorSpec::new("+-", bonds(0.5)).unwrap(),
        OperatorSpec::new("-+", bonds(0.5)).unwrap(),
    ]
}

fn translate(state: u128, sites: usize) -> u128 {
    ((state << 1) & ((1_u128 << sites) - 1)) | (state >> (sites - 1))
}

fn reflect(state: u128, sites: usize) -> u128 {
    let mut result = 0_u128;
    for site in 0..sites {
        result |= ((state >> site) & 1) << (sites - 1 - site);
    }
    result
}

#[test]
fn sector_enumeration_tracks_combinatorial_dimension_across_sizes() {
    for sites in [32, 64, 96] {
        let spin = SpinBasis1D::builder(sites).up(2).build().unwrap();
        let fermion = SpinlessFermionBasis1D::builder(sites)
            .particles(2)
            .build()
            .unwrap();
        let expected = sites * (sites - 1) / 2;
        assert_eq!(spin.len(), expected);
        assert_eq!(fermion.len(), expected);
    }
    for sites in [12, 20, 28] {
        let boson = BosonBasis1D::builder(sites, 3)
            .particles(2)
            .build()
            .unwrap();
        assert_eq!(boson.len(), sites * (sites + 1) / 2);
    }
}

#[test]
fn two_map_general_sectors_remain_sparse_and_invariant() {
    for sites in [6, 8, 10] {
        let translation = ClosureSymmetryMap::new(sites, move |state| {
            Ok((translate(state, sites), Complex64::new(1.0, 0.0)))
        })
        .unwrap();
        let reflection = ClosureSymmetryMap::new(2, move |state| {
            Ok((reflect(state, sites), Complex64::new(1.0, 0.0)))
        })
        .unwrap();
        let basis = GeneralBasis::new(
            SpinBasis1D::builder(sites).up(sites / 2).build().unwrap(),
            SymmetryReducer::new()
                .with_map(translation, 0)
                .with_map(reflection, 0),
        )
        .unwrap();
        let parent_operator = OperatorBuilder::on(basis.parent())
            .terms(periodic_heisenberg(sites))
            .build(MatrixFormat::Csc)
            .unwrap();
        let projector = BasisProjector::from_general(&basis).unwrap();
        assert!(projector
            .preserves_operator_symmetry(&parent_operator, 2.0e-11)
            .unwrap());
        let reduced = OperatorBuilder::on(&basis)
            .terms(periodic_heisenberg(sites))
            .build(MatrixFormat::Csc)
            .unwrap();
        assert!(reduced.nnz() <= basis.len().saturating_mul(2 * sites + 1));
    }
}

#[test]
fn cross_sector_action_and_sparse_algebra_avoid_square_parent_storage() {
    for sites in [10, 14, 18] {
        let source = SpinlessFermionBasis1D::builder(sites)
            .particles(2)
            .build()
            .unwrap();
        let target = SpinlessFermionBasis1D::builder(sites)
            .particles(3)
            .build()
            .unwrap();
        let terms =
            vec![
                OperatorSpec::new("+", (0..sites).map(|site| Coupling::new(1.0, vec![site])))
                    .unwrap(),
            ];
        let source_dimension = f64::from(u32::try_from(source.len()).unwrap());
        let input = vec![Complex64::new(1.0 / source_dimension.sqrt(), 0.0); source.len()];
        let mut streamed = vec![Complex64::new(0.0, 0.0); target.len()];
        apply_sector_shift(&source, &target, &terms, &input, &mut streamed).unwrap();
        let stored = OperatorBuilder::between(&source, &target)
            .terms(terms)
            .build(MatrixFormat::Csc)
            .unwrap();
        let mut expected = vec![Complex64::new(0.0, 0.0); target.len()];
        stored.apply(&input, &mut expected).unwrap();
        assert_eq!(streamed, expected);
        assert!(stored.memory_bytes() < source.len() * target.len() * 16);
    }

    let dimension = 32_768;
    let diagonal = Operator::from_triplets(
        dimension,
        dimension,
        (0..dimension).map(|index| (index, index, Complex64::new(2.0, 0.0))),
        MatrixFormat::Csc,
    )
    .unwrap();
    let polynomial = diagonal.pow(3).unwrap().add(&diagonal).unwrap();
    assert_eq!(polynomial.nnz(), dimension);
    assert!(polynomial.memory_bytes() < 2_000_000);
}

#[test]
fn wide_integer_and_arbitrary_subsystem_oracles_cover_scale_boundaries() {
    let value = (BigUint::from(1_u8) << 700) + (BigUint::from(1_u8) << 193) + BigUint::from(11_u8);
    let encoded: U1024 = state_from_biguint(&value).unwrap();
    assert_eq!(state_to_biguint(encoded), value);
    let wide = WideSpinBasis1024::new(700, Some(1), false).unwrap();
    assert_eq!(wide.len(), 700);

    let amplitude = 1.0 / 2.0_f64.sqrt();
    let mut ghz = vec![Complex64::new(0.0, 0.0); 16];
    ghz[0] = Complex64::new(amplitude, 0.0);
    ghz[15] = Complex64::new(amplitude, 0.0);
    let reduced = partial_trace_subsystem(&ghz, &[2, 2, 2, 2], &[0, 2]).unwrap();
    assert!((reduced[0].re - 0.5).abs() < 1.0e-14);
    assert!((reduced[15].re - 0.5).abs() < 1.0e-14);
    let density: Vec<_> = ghz
        .iter()
        .flat_map(|left| ghz.iter().map(move |right| *left * right.conj()))
        .collect();
    let entropy = entanglement_entropy_density_subsystem(
        &density,
        &[2, 2, 2, 2],
        &[0, 2],
        EntropyOrder::VonNeumann,
    )
    .unwrap();
    assert!((entropy - 2.0_f64.ln()).abs() < 2.0e-13);
}
