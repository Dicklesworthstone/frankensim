//! E5.3b-i calibration sweep (bead wf-root-guzez.6.5.1): drive the FULL
//! nonlinear Dec-17 lifecycle with candidate pilot gains and report the
//! flight-class metrics as JSONL — the evidence a family REGISTRATION
//! cites. Deterministic (fixed seed per candidate); release-mode driver.
//! Repro: cargo run -p fs-flyer --release --bin pilot_nonlinear_sweep

use fs_flyer::pilot::{PilotGains, pilot_family_v1};
use fs_flyer::simloop::{Phase, PilotMode, SimLoop, TerminalEvent, dec17_scenario};

struct FlightClass {
    terminal: &'static str,
    liftoff_tick: Option<u64>,
    end_tick: u64,
    airborne_s: f64,
    undulations: u32,
    max_abs_theta: f64,
    x_end_m: f64,
}

fn classify(mut sim: SimLoop) -> FlightClass {
    let mut liftoff: Option<u64> = None;
    let mut undulations = 0u32;
    let mut last_q_sign = 0i8;
    let mut max_abs_theta: f64 = 0.0;
    let mut x_end = 0.0;
    let mut end_tick = 0;
    let terminal;
    loop {
        match sim.step(None) {
            Ok(out) => {
                x_end = out.x_m;
                end_tick = out.tick;
                if matches!(out.phase, Phase::Airborne) {
                    if liftoff.is_none() {
                        liftoff = Some(out.tick);
                    }
                    max_abs_theta = max_abs_theta.max(out.theta_rad.abs());
                    // Undulation = pitch-rate sign flips while airborne
                    // (two flips per full porpoise cycle; report flips/2).
                    let s = if out.q_rad_s > 1e-3 {
                        1i8
                    } else if out.q_rad_s < -1e-3 {
                        -1i8
                    } else {
                        0
                    };
                    if s != 0 && last_q_sign != 0 && s != last_q_sign {
                        undulations += 1;
                    }
                    if s != 0 {
                        last_q_sign = s;
                    }
                }
                if let Phase::Ended(e) = out.phase {
                    terminal = match e {
                        TerminalEvent::GroundContact => "ground-contact",
                        TerminalEvent::RailEndWithoutLift => "rail-end-without-lift",
                        TerminalEvent::MaxTicks => "max-ticks",
                        TerminalEvent::EnvelopeExceeded => "envelope-exceeded",
                    };
                    break;
                }
            }
            Err(e) => {
                println!("{{\"suite\":\"wf-pilot-sweep\",\"error\":\"{}\"}}", e.code);
                terminal = "hard-refusal";
                break;
            }
        }
    }
    let airborne_s = liftoff.map_or(0.0, |l| (end_tick.saturating_sub(l)) as f64 / 120.0);
    FlightClass {
        terminal,
        liftoff_tick: liftoff,
        end_tick,
        airborne_s,
        undulations: undulations / 2,
        max_abs_theta,
        x_end_m: x_end,
    }
}

fn main() {
    // Phase 2 (registration robustness): the REGISTERED member 3 across
    // seeds — its own remnant stream (tile 3), the real Historical path.
    let _ = pilot_family_v1(3).expect("member 3 registered");
    for seed in 1900u64..1916 {
        let spec = dec17_scenario(seed, PilotMode::Historical(3));
        let sim = SimLoop::init(spec).expect("init");
        let c = classify(sim);
        println!(
            "{{\"suite\":\"wf-pilot-sweep\",\"member\":3,\"seed\":{seed},\"terminal\":\"{}\",\"liftoff_tick\":{},\"end_tick\":{},\"airborne_s\":{:.2},\"undulations\":{},\"max_abs_theta\":{:.3},\"x_end_m\":{:.1}}}",
            c.terminal,
            c.liftoff_tick.map_or(-1i64, |t| t as i64),
            c.end_tick,
            c.airborne_s,
            c.undulations,
            c.max_abs_theta,
            c.x_end_m,
        );
    }
}
