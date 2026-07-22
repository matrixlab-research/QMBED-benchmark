use std::error::Error;
use std::fmt::{Display, Formatter};

use quspin_rust_verify::{LinearOperator, MatrixFormat};

#[derive(Debug)]
struct DimensionError;

impl Display for DimensionError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("dimension mismatch")
    }
}
impl Error for DimensionError {}

struct ScaleOperator {
    shape: (usize, usize),
    factor: f64,
}

impl LinearOperator for ScaleOperator {
    type Scalar = f64;
    type Error = DimensionError;

    fn shape(&self) -> (usize, usize) {
        self.shape
    }

    fn format(&self) -> MatrixFormat {
        MatrixFormat::MatrixFree
    }

    fn apply(
        &self,
        input: &[Self::Scalar],
        output: &mut [Self::Scalar],
    ) -> Result<(), Self::Error> {
        if input.len() != self.shape.1 || output.len() != self.shape.0 {
            return Err(DimensionError);
        }
        for (target, source) in output.iter_mut().zip(input) {
            *target = self.factor * source;
        }
        Ok(())
    }
}

#[test]
fn one_linear_operator_contract_covers_stored_and_matrix_free_algorithms() {
    let operator = ScaleOperator {
        shape: (3, 3),
        factor: 2.0,
    };
    let mut output = [0.0; 3];
    operator.apply(&[1.0, -2.0, 0.5], &mut output).unwrap();
    assert_eq!(operator.shape(), (3, 3));
    assert_eq!(operator.format(), MatrixFormat::MatrixFree);
    assert!(output
        .iter()
        .zip([2.0, -4.0, 1.0])
        .all(|(actual, expected)| (*actual - expected).abs() < f64::EPSILON));
}

#[test]
fn linear_operator_contract_supports_rectangular_sector_changes() {
    let operator = ScaleOperator {
        shape: (2, 3),
        factor: -0.5,
    };
    let mut output = [0.0; 2];
    operator.apply(&[2.0, -4.0, 8.0], &mut output).unwrap();
    assert_eq!(operator.shape(), (2, 3));
    assert!(output
        .iter()
        .zip([-1.0, 2.0])
        .all(|(actual, expected)| (*actual - expected).abs() < f64::EPSILON));
}
