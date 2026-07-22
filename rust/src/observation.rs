use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{Display, Formatter};

use crate::MatrixFormat;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Invariant {
    Finite(&'static str),
    AtMost(&'static str, f64),
    AtLeast(&'static str, f64),
    Between(&'static str, f64, f64),
    Equal(&'static str, f64),
}

impl Invariant {
    #[must_use]
    pub const fn metric(self) -> &'static str {
        match self {
            Self::Finite(metric)
            | Self::AtMost(metric, _)
            | Self::AtLeast(metric, _)
            | Self::Between(metric, _, _)
            | Self::Equal(metric, _) => metric,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Observation {
    pub storage: MatrixFormat,
    pub metrics: BTreeMap<String, f64>,
    pub fingerprint: Option<String>,
}

impl Observation {
    #[must_use]
    pub fn new(storage: MatrixFormat) -> Self {
        Self {
            storage,
            metrics: BTreeMap::new(),
            fingerprint: None,
        }
    }

    #[must_use]
    pub fn metric(mut self, name: impl Into<String>, value: f64) -> Self {
        self.metrics.insert(name.into(), value);
        self
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ValidationError {
    message: String,
}

impl ValidationError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl Display for ValidationError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl Error for ValidationError {}

/// Checks all named scientific invariants supplied by a workflow case.
///
/// # Errors
///
/// Returns [`ValidationError`] if a metric is missing, non-finite, or outside
/// the required bound.
pub fn validate_observation(
    case_id: &str,
    observation: &Observation,
    invariants: &[Invariant],
) -> Result<(), ValidationError> {
    for invariant in invariants {
        let metric = invariant.metric();
        let value = observation.metrics.get(metric).ok_or_else(|| {
            ValidationError::new(format!(
                "{case_id}: observation omitted required metric `{metric}`"
            ))
        })?;
        if !value.is_finite() {
            return Err(ValidationError::new(format!(
                "{case_id}: metric `{metric}` is not finite"
            )));
        }
        let valid = match *invariant {
            Invariant::Finite(_) => true,
            Invariant::AtMost(_, upper) => *value <= upper,
            Invariant::AtLeast(_, lower) => *value >= lower,
            Invariant::Between(_, lower, upper) => *value >= lower && *value <= upper,
            Invariant::Equal(_, expected) => {
                (*value - expected).abs() <= f64::EPSILON * expected.abs().max(1.0)
            }
        };
        if !valid {
            return Err(ValidationError::new(format!(
                "{case_id}: metric `{metric}`={value} violates {invariant:?}"
            )));
        }
    }
    Ok(())
}
