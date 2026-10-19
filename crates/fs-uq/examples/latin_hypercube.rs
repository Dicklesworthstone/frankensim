//! A bounded steady thermal-resistance study with independent uncertain inputs.
//! Run: cargo run -p fs-uq --example latin_hypercube
use fs_uq::{CorrelationModel, LhsExecution, LhsPlan, ParameterUncertainty, UqStatus};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let plan = LhsPlan {
        target_qoi: "steady-body-temperature".into(),
        compliance_threshold: Some(330.0),
        parameters: vec![
            ParameterUncertainty::uniform("ambient", 293.0, 303.0, "K"),
            ParameterUncertainty::uniform("heat-load", 40.0, 60.0, "W"),
            ParameterUncertainty::uniform("conductance", 1.5, 2.5, "W/K"),
        ],
        correlation: CorrelationModel::Independent,
        replicates: 8,
        samples_per_replicate: 127,
        seed: 0x4c48_5301,
    };
    let mut execution = LhsExecution::new(&plan)?;
    let report = execution.advance(usize::MAX, || false, |parameters| {
        Ok::<_, &'static str>(parameters[0] + parameters[1] / parameters[2])
    });
    if report.status != UqStatus::Complete {
        return Err(format!("study refused: {:?}", report.rejection_reason).into());
    }
    let estimate = report.estimate.as_ref().ok_or("no complete replicate estimate")?;
    println!("{}: mean={} K, between-replicate standard_error={:?} K",
        report.qoi_name, estimate.mean, estimate.standard_error);
    println!("complete_replicates={} accepted_samples={} empirical_compliance={:?}",
        report.completed_replicates, report.samples_accepted, report.compliance);
    println!("Estimated thermal-resistance model only; no model-error or optional-stopping bound.");
    Ok(())
}
