use quspin_rust_verify::candidate::QuSpinAdapter;
use quspin_rust_verify::{paper_workflows, validate_observation, WorkflowBackend};

#[test]
#[ignore = "paper-scale workflows; exercised in release mode"]
fn pinned_candidate_passes_all_twelve_paper_workflows() {
    let mut backend = QuSpinAdapter;
    for case in paper_workflows() {
        let observation = backend
            .run(&case)
            .unwrap_or_else(|error| panic!("{} backend failed: {error}", case.case_id));
        validate_observation(case.case_id, &observation, case.invariants)
            .unwrap_or_else(|error| panic!("{} validation failed: {error}", case.case_id));
    }
}
