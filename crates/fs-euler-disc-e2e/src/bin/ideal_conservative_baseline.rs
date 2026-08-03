//! Deterministic executable for the ideal no-slip Euler-disc baseline rung.

use fs_euler_disc_e2e::{SquatDiscInput, run_ideal_conservative_baseline};

fn main() {
    let output = run_ideal_conservative_baseline(SquatDiscInput::nominal());
    println!("{}", output.structured_output());
}
