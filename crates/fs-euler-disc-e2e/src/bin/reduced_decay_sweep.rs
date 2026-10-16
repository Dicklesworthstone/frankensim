//! Deterministic numerical-reference sweep for the public reduced-decay API.

use fs_euler_disc_e2e::{
    ReducedDecayInput, channel_crossover_diagnostic, refinement_evidence, run_reduced_decay,
    structured_runner_output,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = ReducedDecayInput::nominal_reference()?;
    let run = run_reduced_decay(&input)?;
    let refinement = refinement_evidence(&input)?;
    println!("{}", structured_runner_output(&run, &refinement)?);
    println!(
        "schema=reduced-decay-v1 channel_crossover={:?} evidence_scope=encoded-model-only",
        channel_crossover_diagnostic(&input)?,
    );
    Ok(())
}
