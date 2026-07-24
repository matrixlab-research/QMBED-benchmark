use std::collections::HashSet;

use qmbed_verify::{paper_workflows, Capability, WorkflowOrigin};

#[test]
fn catalog_preserves_the_twelve_python_julia_workflows() {
    let cases = paper_workflows();
    assert_eq!(cases.len(), 12);
    assert_eq!(
        cases
            .iter()
            .filter(|case| case.origin == WorkflowOrigin::OriginalBaseline)
            .count(),
        6
    );
    assert_eq!(
        cases
            .iter()
            .filter(|case| case.origin == WorkflowOrigin::LkmExtension)
            .count(),
        6
    );

    let ids = cases
        .iter()
        .map(|case| case.case_id)
        .collect::<HashSet<_>>();
    assert_eq!(ids.len(), cases.len());
    assert!(ids.contains("paper_mbl_shift_invert_l14"));
    assert!(ids.contains("paper_pxp_revival_l24"));
    assert!(ids.contains("paper_particle_addition_6x3"));
}

#[test]
fn catalog_forces_basis_and_algorithm_breadth() {
    let cases = paper_workflows();
    let capabilities = cases
        .iter()
        .flat_map(|case| case.capabilities.iter().copied())
        .collect::<HashSet<_>>();

    for required in [
        Capability::SpinBasis,
        Capability::BosonBasis,
        Capability::SpinlessFermionBasis,
        Capability::SpinfulFermionBasis,
        Capability::UserBasis,
        Capability::SymmetrySector,
        Capability::CscAssembly,
        Capability::ShiftInvert,
        Capability::KrylovEvolution,
        Capability::SubspaceTracking,
        Capability::KrylovSpectrum,
        Capability::CrossSectorOperator,
    ] {
        assert!(capabilities.contains(&required), "missing {required:?}");
    }
}

#[test]
fn every_case_has_executable_invariants() {
    for case in paper_workflows() {
        assert!(
            !case.invariants.is_empty(),
            "{} has no invariant",
            case.case_id
        );
        assert!(!case.parameters.is_empty());
        assert!(!case.capabilities.is_empty());
    }
}
