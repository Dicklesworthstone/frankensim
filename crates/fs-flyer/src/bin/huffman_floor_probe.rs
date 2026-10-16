//! E5.4 scratch probe (bead wf-root-guzez.6.7): find the LOWEST
//! headwind at which the certified coupled aero core carries a Huffman
//! catapult rail run to liftoff. JSONL per candidate; release-mode.
//! Repro: cargo run -p fs-flyer --release --bin huffman_floor_probe

use fs_flyer::simloop::{Phase, PilotMode, SimLoop, huffman_scenario};

fn main() {
    for headwind in [2.0f64, 3.0, 4.0, 5.0, 6.0, 8.0] {
        let mut spec = huffman_scenario(1904, PilotMode::FixedControls);
        spec.headwind_mps = headwind;
        let mut sim = SimLoop::init(spec).expect("init");
        let mut outcome = "running";
        let mut tick = 0;
        let mut x = 0.0;
        loop {
            match sim.step(None) {
                Ok(out) => {
                    tick = out.tick;
                    x = out.x_m;
                    if matches!(out.phase, Phase::Airborne) {
                        outcome = "liftoff";
                        break;
                    }
                    if let Phase::Ended(e) = out.phase {
                        outcome = match e {
                            fs_flyer::simloop::TerminalEvent::EnvelopeExceeded => "envelope",
                            fs_flyer::simloop::TerminalEvent::RailEndWithoutLift => "rail-end",
                            _ => "other",
                        };
                        break;
                    }
                }
                Err(e) => {
                    println!("{{\"headwind\":{headwind},\"error\":\"{}\"}}", e.code);
                    outcome = "refusal";
                    break;
                }
            }
        }
        println!(
            "{{\"suite\":\"wf-huffman-probe\",\"headwind\":{headwind},\"outcome\":\"{outcome}\",\"tick\":{tick},\"x_m\":{x:.2}}}"
        );
    }
}
