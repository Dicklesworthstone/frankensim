//! E6.2 bisection probe. Repro: cargo run -p fs-flyer --release --bin tick0probe
use fs_flyer::simloop::{PilotMode, SimLoop, dec17_scenario};
fn main() {
    let mut spec = dec17_scenario(1903, PilotMode::FixedControls);
    spec.max_ticks = 120;
    let sim = SimLoop::init(spec).expect("init");
    let t = &sim.tick0().trim;
    println!(
        "{{\"lane\":\"native\",\"v\":\"{}\",\"alpha\":\"{}\",\"dc\":\"{}\",\"omega\":\"{}\",\"iters\":{}}}",
        t.v_mps.to_bits(),
        t.alpha_rad.to_bits(),
        t.delta_canard_rad.to_bits(),
        t.omega_prop_rad_s.to_bits(),
        t.iterations
    );
}
