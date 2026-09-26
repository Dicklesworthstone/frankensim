//! A bounded, ideal point-mass bounce; not a calibrated contact model.
//! Run: cargo run -p fs-time --example rk45_impact
use fs_time::{
    AdaptiveState, EventAdvance, EventDirection, EventOptions, InitialEvent,
    PiController, rk45_until_event,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // SI units: seconds, metres, metres/second. Acceleration is explicit.
    let gravity = 9.81;
    let restitution = 0.8;
    let rhs = |_: f64, u: &[f64], out: &mut [f64]| {
        out[0] = u[1];
        out[1] = -gravity;
    };
    let guard = |_: f64, u: &[f64]| u[0];
    let pi = PiController::default();
    let mut state = AdaptiveState::new(0.0, &[10.0, 0.0], 0.1);
    let mut options = EventOptions {
        direction: EventDirection::Falling,
        initial: InitialEvent::Report,
        max_step: 0.25,
        scan_substeps: 8,
        time_tolerance: 1e-10,
        max_iterations: 80,
    };
    for impact in 1..=3 {
        match rk45_until_event(&mut state, &rhs, &guard, 10.0, 1e-9, 1e-11,
                              &pi, 10_000, &options, &mut || false)? {
            EventAdvance::Event(event) => {
                println!("impact={impact} time_s={:.10} velocity_m_s={:.10} bracket_s={:?}",
                         event.time, state.u[1], event.bracket);
                // The caller, not the root locator, owns the reset law.
                state.u[0] = 0.0;
                state.u[1] *= -restitution;
                options.initial = InitialEvent::Ignore;
            }
            EventAdvance::Stopped(report) => {
                return Err(format!("integration stopped before impact {impact}: {report:?}").into());
            }
        }
    }
    Ok(())
}
