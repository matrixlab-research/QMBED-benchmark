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
    dimension: usize,
    factor: f64,
}

impl LinearOperator for ScaleOperator {
    type Scalar = f64;
    type Error = DimensionError;

    fn dimension(&self) -> usize {
        self.dimension
    }

    fn format(&self) -> MatrixFormat {
        MatrixFormat::MatrixFree
    }

    fn apply(
        &self,
        input: &[Self::Scalar],
        output: &mut [Self::Scalar],
    ) -> Result<(), Self::Error> {
        if input.len() != self.dimension || output.len() != self.dimension {
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
        dimension: 3,
        factor: 2.0,
    };
    let mut output = [0.0; 3];
    operator.apply(&[1.0, -2.0, 0.5], &mut output).unwrap();
    assert_eq!(operator.dimension(), 3);
    assert_eq!(operator.format(), MatrixFormat::MatrixFree);
    assert!(output
        .iter()
        .zip([2.0, -4.0, 1.0])
        .all(|(actual, expected)| (*actual - expected).abs() < f64::EPSILON));
}
