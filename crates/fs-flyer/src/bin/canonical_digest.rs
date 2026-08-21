//! E6.2 canonical-scenario digest probe (bead wf-root-guzez.7.2):
//! print the 1-second canonical digest (dec17 fixed, seed 1903, 120
//! ticks) for lane comparison. Repro:
//! cargo run -p fs-flyer --release --bin canonical_digest

use fs_flyer::simloop::{Phase, PilotMode, SimLoop, dec17_scenario};

fn main() {
    let mut spec = dec17_scenario(1903, PilotMode::FixedControls);
    spec.max_ticks = 120;
    let mut sim = SimLoop::init(spec).expect("init");
    loop {
        let out = sim.step(None).expect("step");
        if let Phase::Ended(_) = out.phase {
            break;
        }
    }
    println!(
        "{{\"suite\":\"wf-canonical\",\"lane\":\"native\",\"digest\":\"{}\"}}",
        sim.digest_hex()
    );
}
