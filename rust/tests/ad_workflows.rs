use std::collections::BTreeSet;

use qmbed_benchmark::ad_workflows::{ad_workflows, benchmark_ad_workflow, AdScale};

#[test]
fn catalog_contains_twelve_distinct_parameterized_workflows() {
    let workflows = ad_workflows(AdScale::Ci).unwrap();
    assert_eq!(workflows.len(), 12);
    assert_eq!(
        workflows
            .iter()
            .map(|workflow| workflow.case_id)
            .collect::<BTreeSet<_>>()
            .len(),
        12
    );
    assert!(workflows
        .iter()
        .all(|workflow| workflow.parameters.values().len() >= 2));
}

#[test]
#[ignore = "release-mode AD acceptance; exercised explicitly in CI"]
fn ci_workflows_match_finite_differences_and_report_algorithmic_advantage() {
    for workflow in ad_workflows(AdScale::Ci).unwrap() {
        let row = benchmark_ad_workflow(&workflow, 0, 1)
            .unwrap_or_else(|error| panic!("{} failed: {error}", workflow.case_id));
        assert!(row.maximum_absolute_error <= 2.0e-5 || row.maximum_relative_error <= 2.0e-4);
        assert!(row.spectral_gap > 0.0);
        assert!(row.residual <= 1.0e-8);
        assert_eq!(row.analytic_eigensolves, 1);
        assert_eq!(row.finite_difference_eigensolves, 2 * row.parameters);
    }
}

#[test]
#[ignore = "paper-scale AD workflows; exercised by the scheduled benchmark job"]
fn paper_scale_workflows_match_finite_differences() {
    for workflow in ad_workflows(AdScale::Paper).unwrap() {
        benchmark_ad_workflow(&workflow, 0, 1)
            .unwrap_or_else(|error| panic!("{} failed: {error}", workflow.case_id));
    }
}
