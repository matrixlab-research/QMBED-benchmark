use std::error::Error;
use std::fmt::{Display, Formatter};
use std::time::Instant;

use crate::{validate_observation, Observation, ValidationError, WorkflowCase};

#[allow(clippy::missing_errors_doc)]
pub trait WorkflowBackend {
    type Error: Error + Send + Sync + 'static;

    fn language(&self) -> &'static str;
    fn run(&mut self, case: &WorkflowCase) -> Result<Observation, Self::Error>;
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BenchmarkOptions {
    pub warmups: usize,
    pub samples: usize,
}

impl Default for BenchmarkOptions {
    fn default() -> Self {
        Self {
            warmups: 1,
            samples: 5,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BenchmarkRow {
    pub suite: &'static str,
    pub case_id: &'static str,
    pub family_id: &'static str,
    pub benchmark: &'static str,
    pub parameters: &'static str,
    pub language: &'static str,
    pub comparison: &'static str,
    pub storage: &'static str,
    pub validation: &'static str,
    pub samples: usize,
    pub iterations_per_sample: usize,
    pub minimum_seconds: f64,
    pub p05_seconds: f64,
    pub p25_seconds: f64,
    pub median_seconds: f64,
    pub mean_seconds: f64,
    pub p75_seconds: f64,
    pub p95_seconds: f64,
    pub maximum_seconds: f64,
    pub stddev_seconds: f64,
    pub raw_samples_seconds: Vec<f64>,
}

impl BenchmarkRow {
    pub const CSV_HEADER: &'static str = "suite,case_id,family_id,benchmark,parameters,language,comparison,storage,validation,samples,iterations_per_sample,min_seconds,p05_seconds,p25_seconds,median_seconds,mean_seconds,p75_seconds,p95_seconds,max_seconds,stddev_seconds,raw_samples_seconds";

    #[must_use]
    pub fn to_csv_record(&self) -> String {
        let raw = self
            .raw_samples_seconds
            .iter()
            .map(|value| format!("{value:.12}"))
            .collect::<Vec<_>>()
            .join(";");
        [
            self.suite.to_string(),
            self.case_id.to_string(),
            self.family_id.to_string(),
            csv_escape(self.benchmark),
            csv_escape(self.parameters),
            self.language.to_string(),
            self.comparison.to_string(),
            self.storage.to_string(),
            self.validation.to_string(),
            self.samples.to_string(),
            self.iterations_per_sample.to_string(),
            format!("{:.12}", self.minimum_seconds),
            format!("{:.12}", self.p05_seconds),
            format!("{:.12}", self.p25_seconds),
            format!("{:.12}", self.median_seconds),
            format!("{:.12}", self.mean_seconds),
            format!("{:.12}", self.p75_seconds),
            format!("{:.12}", self.p95_seconds),
            format!("{:.12}", self.maximum_seconds),
            format!("{:.12}", self.stddev_seconds),
            raw,
        ]
        .join(",")
    }
}

#[derive(Debug)]
pub enum BenchmarkError<E> {
    Backend(E),
    Validation(ValidationError),
    InvalidOptions(&'static str),
}

impl<E: Display> Display for BenchmarkError<E> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Backend(error) => write!(formatter, "candidate backend failed: {error}"),
            Self::Validation(error) => write!(formatter, "workflow validation failed: {error}"),
            Self::InvalidOptions(message) => formatter.write_str(message),
        }
    }
}

impl<E: Error + 'static> Error for BenchmarkError<E> {}

/// Runs correctness preflights before collecting warm workflow samples.
///
/// # Errors
///
/// Returns an error when the backend fails, an observation violates a physical
/// invariant, or the sampling configuration is empty.
pub fn benchmark_suite<B: WorkflowBackend>(
    backend: &mut B,
    cases: &[WorkflowCase],
    options: BenchmarkOptions,
) -> Result<Vec<BenchmarkRow>, BenchmarkError<B::Error>> {
    if options.samples == 0 {
        return Err(BenchmarkError::InvalidOptions(
            "benchmark samples must be greater than zero",
        ));
    }
    let mut rows = Vec::with_capacity(cases.len());
    for case in cases {
        for _ in 0..options.warmups {
            let observation = backend.run(case).map_err(BenchmarkError::Backend)?;
            validate_observation(case.case_id, &observation, case.invariants)
                .map_err(BenchmarkError::Validation)?;
        }

        let mut samples = Vec::with_capacity(options.samples);
        let mut storage = "unspecified";
        for _ in 0..options.samples {
            let started = Instant::now();
            let observation = backend.run(case).map_err(BenchmarkError::Backend)?;
            let elapsed = started.elapsed().as_secs_f64();
            validate_observation(case.case_id, &observation, case.invariants)
                .map_err(BenchmarkError::Validation)?;
            storage = observation.storage.as_str();
            samples.push(elapsed);
        }
        rows.push(summarize(case, backend.language(), storage, samples));
    }
    Ok(rows)
}

#[allow(clippy::cast_precision_loss)]
fn summarize(
    case: &WorkflowCase,
    language: &'static str,
    storage: &'static str,
    mut samples: Vec<f64>,
) -> BenchmarkRow {
    samples.sort_by(f64::total_cmp);
    let count = samples.len() as f64;
    let mean = samples.iter().sum::<f64>() / count;
    let variance = samples
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / count;
    BenchmarkRow {
        suite: "paper",
        case_id: case.case_id,
        family_id: case.family_id,
        benchmark: case.name,
        parameters: case.parameters,
        language,
        comparison: "current_backend",
        storage,
        validation: "physical_invariant",
        samples: samples.len(),
        iterations_per_sample: 1,
        minimum_seconds: samples[0],
        p05_seconds: percentile(&samples, 0.05),
        p25_seconds: percentile(&samples, 0.25),
        median_seconds: percentile(&samples, 0.50),
        mean_seconds: mean,
        p75_seconds: percentile(&samples, 0.75),
        p95_seconds: percentile(&samples, 0.95),
        maximum_seconds: samples[samples.len() - 1],
        stddev_seconds: variance.sqrt(),
        raw_samples_seconds: samples,
    }
}

// Benchmark sample counts are deliberately tiny (five for the paper suite),
// so these bounded index conversions cannot lose meaningful precision.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]
fn percentile(sorted: &[f64], probability: f64) -> f64 {
    let position = probability * (sorted.len() - 1) as f64;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    if lower == upper {
        return sorted[lower];
    }
    let weight = position - lower as f64;
    sorted[lower] * (1.0 - weight) + sorted[upper] * weight
}

fn csv_escape(value: &str) -> String {
    if value.contains([',', '"', '\n']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}
