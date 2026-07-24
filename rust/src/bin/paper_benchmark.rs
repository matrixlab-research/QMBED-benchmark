use std::env;
use std::error::Error;

use qmbed_verify::candidate::QmbedAdapter;
use qmbed_verify::{benchmark_suite, paper_workflows, BenchmarkOptions, BenchmarkRow};

fn count_from_environment(name: &str, default: usize) -> Result<usize, Box<dyn Error>> {
    match env::var(name) {
        Ok(value) => Ok(value.parse()?),
        Err(env::VarError::NotPresent) => Ok(default),
        Err(error) => Err(Box::new(error)),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let options = BenchmarkOptions {
        warmups: count_from_environment("QMBED_RUST_WARMUPS", 1)?,
        samples: count_from_environment("QMBED_RUST_SAMPLES", 5)?,
    };
    let rows = benchmark_suite(&mut QmbedAdapter, &paper_workflows(), options)?;
    println!("{}", BenchmarkRow::CSV_HEADER);
    for row in rows {
        println!("{}", row.to_csv_record());
    }
    Ok(())
}
