use std::convert::Infallible;

use qmbed_benchmark::{
    benchmark_suite, paper_workflows, BenchmarkOptions, Invariant, MatrixFormat, Observation,
    WorkflowBackend,
};

struct ContractBackend;

impl WorkflowBackend for ContractBackend {
    type Error = Infallible;

    fn language(&self) -> &'static str {
        "rust"
    }

    fn run(&mut self, case: &qmbed_benchmark::WorkflowCase) -> Result<Observation, Self::Error> {
        let mut observation = Observation::new(MatrixFormat::Csc);
        for invariant in case.invariants {
            let value = match *invariant {
                Invariant::Finite(_) => 0.0,
                Invariant::AtMost(_, upper) => upper * 0.5,
                Invariant::AtLeast(_, lower) => {
                    if lower >= 0.0 {
                        lower + 1.0
                    } else {
                        0.0
                    }
                }
                Invariant::Between(_, lower, upper) => (lower + upper) * 0.5,
                Invariant::Equal(_, expected) => expected,
            };
            observation = observation.metric(invariant.metric(), value);
        }
        Ok(observation)
    }
}
#[test]
fn rust_rows_match_the_existing_paper_csv_contract() {
    let rows = benchmark_suite(
        &mut ContractBackend,
        &paper_workflows(),
        BenchmarkOptions {
            warmups: 1,
            samples: 2,
        },
    )
    .expect("contract backend should satisfy all invariants");

    assert_eq!(rows.len(), 12);
    assert!(rows.iter().all(|row| row.language == "rust"));
    assert!(rows.iter().all(|row| row.suite == "paper"));
    assert!(rows.iter().all(|row| row.samples == 2));

    let csv = format!(
        "{}\n{}\n",
        qmbed_benchmark::BenchmarkRow::CSV_HEADER,
        rows[0].to_csv_record()
    );
    assert!(csv
        .lines()
        .nth(1)
        .is_some_and(|row| row.starts_with("rust,")));
    assert!(csv.contains("raw_samples_seconds"));
    assert!(csv.contains("stdev_seconds"));
    assert!(!csv.contains("stddev_seconds"));
    assert!(csv.contains(",workflow,end_to_end,csc,true,"));
}

#[test]
fn invalid_observations_fail_before_timing_is_reported() {
    struct InvalidBackend;

    impl WorkflowBackend for InvalidBackend {
        type Error = Infallible;

        fn language(&self) -> &'static str {
            "rust"
        }

        fn run(
            &mut self,
            _case: &qmbed_benchmark::WorkflowCase,
        ) -> Result<Observation, Self::Error> {
            Ok(Observation::new(MatrixFormat::Csc).metric("residual", f64::NAN))
        }
    }

    let first = &paper_workflows()[0..1];
    let error = benchmark_suite(
        &mut InvalidBackend,
        first,
        BenchmarkOptions {
            warmups: 0,
            samples: 1,
        },
    )
    .expect_err("invalid physical output must reject the timing row");
    assert!(error.to_string().contains("not finite"));
}
