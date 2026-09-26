//! Fit restitution to a synthetic third-impact time using saltation-aware
//! derivatives. Ideal point-mass impacts, not calibrated contact physics.
//! Run: cargo run -p fs-time --example rk45_impact_fit -- 5.54003131687

use fs_time::{EventDirection, InitialEvent, PiController};
use fs_time::adaptive::events::multiple::{EventOrder, EventSetHit, EventSetOptions, EventSpec};
use fs_time::adaptive::events::multiple::run::{
    HybridConfig, HybridReset, HybridStop, HybridSystem, ResetAction,
};
use fs_time::hybrid_sensitivity::{ForwardSensitivity, HybridDerivatives};

type Error = Box<dyn std::error::Error>;
const FLOOR: [EventSpec; 1] = [EventSpec {
    id: 1, direction: EventDirection::Falling, initial: InitialEvent::Report,
}];
struct Impact { restitution: f64 }
impl HybridSystem for Impact {
    type Mode = usize;
    fn rhs(&self, _: &usize, _: f64, u: &[f64], out: &mut [f64]) {
        out[0] = u[1]; out[1] = -9.81;
    }
    fn events(&self, _: &usize) -> &[EventSpec] { &FLOOR }
    fn guard(&self, _: &usize, _: u64, _: f64, u: &[f64]) -> f64 { u[0] }
    fn reset(&self, mode: &usize, hit: &EventSetHit, u: &[f64]) -> Result<HybridReset<usize>, String> {
        Ok(HybridReset { mode: mode + 1, state: vec![0.0, -self.restitution * u[1]],
            action: if *mode == 2 { ResetAction::Terminate } else { ResetAction::Continue },
            consumed_ids: vec![hit.selected.id] })
    }
}
impl HybridDerivatives for Impact {
    fn rhs_tangent(&self, _: &usize, _: f64, _: &[f64], _: usize, s: &[f64], out: &mut [f64]) {
        out[0] = s[1]; out[1] = 0.0;
    }
    fn guard_gradient(&self, _: &usize, _: u64, _: f64, _: &[f64], out: &mut [f64]) -> Result<f64, String> {
        out.copy_from_slice(&[1.0, 0.0]); Ok(0.0)
    }
    fn guard_parameter(&self, _: &usize, _: u64, _: f64, _: &[f64], _: usize) -> Result<f64, String> { Ok(0.0) }
    fn reset_tangent(&self, _: &usize, _: &EventSetHit, u: &[f64], _: &HybridReset<usize>,
        _: usize, s: &[f64], _: f64, out: &mut [f64]) -> Result<(), String>
    { out[0] = 0.0; out[1] = -self.restitution * s[1] - u[1]; Ok(()) }
    fn crossing_speed_floor(&self, _: &usize, _: u64) -> f64 { 1e-8 }
}

fn evaluate(restitution: f64) -> Result<(f64, f64), Error> {
    let model = Impact { restitution };
    let forward = ForwardSensitivity::new(&model, 2, 1, 16)?;
    // Direction 0 varies restitution only. SI: metres, seconds, m/s.
    let mut state = forward.initial(0, 0.0, &[10.0, 0.0], 0.1, &[0.0, 0.0])?;
    let cfg = HybridConfig { t_end: 10.0, rtol: 1e-10, atol: 1e-12,
        pi: PiController::default(), max_attempts: 20_000, max_events: 3,
        events: EventSetOptions { max_step: 0.2, scan_substeps: 8,
            time_tolerance: 1e-11, max_iterations: 80, order: EventOrder::RequireSeparated } };
    let report = forward.run(&mut state, &cfg, &mut || false)?;
    if report.stop != HybridStop::Terminated { return Err(format!("incomplete impact solve: {report:?}").into()); }
    let gradient = state.last_event_time_tangents().ok_or("missing terminal-time derivative")?[0];
    Ok((state.trajectory().integration().t, gradient))
}

fn fit(target: f64) -> Result<(f64, usize), Error> {
    if !target.is_finite() || target <= 0.0 { return Err("target time must be positive and finite".into()); }
    let (mut lo, mut hi) = (0.05, 0.95);
    let (low_time, _) = evaluate(lo)?;
    let (high_time, _) = evaluate(hi)?;
    if target < low_time || target > high_time { return Err("target is outside the admitted restitution bracket".into()); }
    if (target - low_time).abs() < 1e-9 { return Ok((lo, 0)); }
    if (target - high_time).abs() < 1e-9 { return Ok((hi, 0)); }
    let mut e = 0.5;
    for iteration in 1..=16 {
        let (time, derivative) = evaluate(e)?;
        let residual = time - target;
        if residual.abs() < 1e-9 { return Ok((e, iteration)); }
        if residual < 0.0 { lo = e; } else { hi = e; }
        let proposal = e - residual / derivative;
        e = if proposal.is_finite() && proposal > lo && proposal < hi { proposal } else { (lo + hi) / 2.0 };
    }
    Err("restitution fit exhausted its iteration budget".into())
}

fn main() -> Result<(), Error> {
    let target = std::env::args().nth(1).map(|s| s.parse::<f64>()).transpose()?.unwrap_or(5.54003131687);
    let (e, iterations) = fit(target)?;
    let (time, gradient) = evaluate(e)?;
    println!("restitution={e:.10} target_s={target:.10} simulated_s={time:.10} dt_de={gradient:.10} iterations={iterations}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn recovers_synthetic_restitution_through_three_event_time_derivatives() {
        let target = (20.0f64 / 9.81).sqrt() * (1.0 + 2.0 * 0.8 + 2.0 * 0.8 * 0.8);
        let (e, iterations) = fit(target).unwrap();
        assert!((e - 0.8).abs() < 1e-8);
        assert!(iterations <= 8);
    }
    #[test]
    fn invalid_or_unreachable_targets_refuse() {
        assert!(fit(f64::NAN).is_err());
        assert!(fit(-1.0).is_err());
        assert!(fit(0.5).is_err());
        assert!(fit(20.0).is_err());
    }
}
