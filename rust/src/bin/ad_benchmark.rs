use std::env;
use std::error::Error;

use qmbed_benchmark::ad_workflows::{ad_workflows, benchmark_ad_workflow, AdBenchmarkRow, AdScale};

fn count(name: &str, default: usize) -> Result<usize, Box<dyn Error>> {
    Ok(env::var(name).map_or_else(|_| Ok(default), |value| value.parse())?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let scale = match env::var("QMBED_AD_SCALE") {
        Ok(value) if value == "paper" => AdScale::Paper,
        Ok(value) if value == "ci" => AdScale::Ci,
        Err(env::VarError::NotPresent) => AdScale::Ci,
        Ok(value) => return Err(format!("unknown QMBED_AD_SCALE {value:?}").into()),
        Err(error) => return Err(Box::new(error)),
    };
    let warmups = count("QMBED_AD_WARMUPS", 1)?;
    let samples = count("QMBED_AD_SAMPLES", 3)?;
    println!("{}", AdBenchmarkRow::CSV_HEADER);
    for case in ad_workflows(scale)? {
        println!(
            "{}",
            benchmark_ad_workflow(&case, warmups, samples)?.to_csv_record()
        );
    }
    Ok(())
}
