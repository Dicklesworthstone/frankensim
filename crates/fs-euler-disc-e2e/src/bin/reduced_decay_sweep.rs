//! Deterministic numerical-reference sweep for the isolated reduced-decay module.

#[path = "../reduced_decay.rs"]
mod reduced_decay;

use reduced_decay::{
    ReducedDecayInput, refinement_evidence, run_reduced_decay, structured_runner_output,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = ReducedDecayInput::nominal_reference()?;
    let run = run_reduced_decay(&input)?;
    let refinement = refinement_evidence(&input)?;
    println!("{}", structured_runner_output(&run, &refinement)?);
    Ok(())
}
