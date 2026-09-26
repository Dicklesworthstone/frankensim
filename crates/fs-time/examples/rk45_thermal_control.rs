//! Gradient-driven heater control for an ideal lumped thermal model.
//!
//! SI inputs: capacity 100 J/K, conductance 2 W/K, ambient 300 K, horizon
//! 60 s, heater bounds 0..60 W, target 310 K. Eight Bernstein coefficients
//! parameterize a smooth power schedule; their convex hull enforces the bounds.
//! J = 0.5*((T_end-310 K)/10 K)^2 + 0.01*mean((power/60 W)^2).
//! The second term is a quadratic control penalty, NOT physical energy.
//! The whole gradient comes from one discrete reverse sweep, not eight
//! finite-difference trajectories. This is an illustrative model, not validated
//! cooling hardware, a constrained-temperature guarantee or a global optimum.
use fs_time::{AdaptiveState, PiController};
use fs_time::adaptive::adjoint::OdeVjp;
use fs_time::adaptive::adjoint::trajectory::{RecordedRk45, RecordingConfig, RecordingStatus, ReplayBudget};

const CAPACITY: f64 = 100.0;
const CONDUCTANCE: f64 = 2.0;
const AMBIENT: f64 = 300.0;
const HORIZON: f64 = 60.0;
const MAX_POWER: f64 = 60.0;
const TARGET: f64 = 310.0;
const SCALE: f64 = 10.0;
const PENALTY: f64 = 0.01;
type Error = Box<dyn std::error::Error>;

fn basis(time: f64) -> [f64; 8] {
    let s = time / HORIZON;
    let mut b = [1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0];
    for (i, value) in b.iter_mut().enumerate() {
        for _ in 0..i { *value *= s; }
        for _ in i..7 { *value *= 1.0 - s; }
    }
    b
}
struct Heater { controls: [f64; 8] }
impl Heater {
    fn power(&self, b: &[f64; 8]) -> f64 {
        b.iter().zip(self.controls).fold(0.0, |sum, (b,p)| b.mul_add(p,sum))
    }
}
impl OdeVjp for Heater {
    fn dimension(&self) -> usize { 2 }
    fn parameter_count(&self) -> usize { 8 }
    fn rhs(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let power = self.power(&basis(t));
        out[0] = (power - CONDUCTANCE*(x[0]-AMBIENT))/CAPACITY;
        out[1] = (power/MAX_POWER)*(power/MAX_POWER)/HORIZON;
    }
    fn rhs_vjp(&self, t: f64, _: &[f64], seed: &[f64], x: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        let b = basis(t);
        let power = self.power(&b);
        x[0] = -CONDUCTANCE/CAPACITY*seed[0]; x[1] = 0.0;
        let multiplier = seed[0]/CAPACITY + seed[1]*2.0*power/(MAX_POWER*MAX_POWER*HORIZON);
        for i in 0..8 { p[i] = b[i]*multiplier; }
        Ok(())
    }
}
struct Evaluation { cost: f64, temperature: f64, gradient: Vec<f64>, replays: usize, checkpoints: usize }
fn evaluate(controls: [f64; 8], derivatives: bool) -> Result<Evaluation, Error> {
    if controls.iter().any(|v| !v.is_finite() || *v < 0.0 || *v > MAX_POWER) {
        return Err("heater controls must lie within the explicit finite power bounds".into());
    }
    let model = Heater { controls };
    let config = RecordingConfig { end: HORIZON, rtol: 1e-10, atol: 1e-12,
        controller: PiController::default(), max_workspace_components: 1024 };
    let mut trajectory = RecordedRk45::new(&model, AdaptiveState::new(0.0, &[AMBIENT,0.0],0.5), config)?;
    let report = trajectory.advance(10000,10000,&mut||false)?;
    if report.status != RecordingStatus::ReachedEnd { return Err(format!("incomplete solve: {report:?}").into()); }
    let state = &trajectory.state().u;
    let temperature = state[0];
    let error = (temperature-TARGET)/SCALE;
    let cost = 0.5*error*error + PENALTY*state[1];
    let mut result = Evaluation { cost, temperature, gradient: Vec::new(), replays: 0, checkpoints: 0 };
    if derivatives {
        let gradient = trajectory.pullback(&[error/SCALE,PENALTY], &[0.0;8],
            ReplayBudget { checkpoints: 32, replayed_steps: 100_000 }, &mut||false)?;
        result.gradient = gradient.parameters;
        result.replays = gradient.replayed_steps;
        result.checkpoints = gradient.peak_checkpoints;
    }
    Ok(result)
}

/// Safeguarded projected descent, deliberately bounded. Iteration exhaustion
/// is reported as such, never labeled convergence. Each candidate gets a fresh
/// accuracy-controlled primal solve; each gradient freezes that solve's mesh.
fn optimize(iterations: usize) -> Result<([f64;8], Evaluation, &'static str), Error> {
    let mut controls = [5.0;8];
    for _ in 0..iterations {
        let current = evaluate(controls,true)?;
        let projected = controls.iter().zip(&current.gradient).fold(0.0f64, |norm,(p,g)| {
            norm.max((p-(p-g).clamp(0.0,MAX_POWER)).abs())
        });
        if projected < 1e-7 { return Ok((controls,current,"projected-gradient tolerance")); }
        let mut step = 5000.0;
        let mut accepted = None;
        for _ in 0..16 {
            let mut candidate = controls;
            let mut directional = 0.0;
            for i in 0..8 {
                candidate[i] = (controls[i]-step*current.gradient[i]).clamp(0.0,MAX_POWER);
                directional += current.gradient[i]*(candidate[i]-controls[i]);
            }
            if directional < 0.0 {
                let trial = evaluate(candidate,false)?;
                if trial.cost <= current.cost+1e-4*directional { accepted=Some(candidate);break; }
            }
            step *= 0.5;
        }
        controls = accepted.ok_or("projected line search exhausted its trial budget")?;
    }
    Ok((controls,evaluate(controls,true)?,"iteration budget exhausted"))
}
fn main() -> Result<(), Error> {
    let initial = evaluate([5.0;8],false)?;
    let (controls,result,status) = optimize(80)?;
    println!("status={status}");
    println!("initial_cost={:.10} final_cost={:.10} final_temperature_K={:.8}",initial.cost,result.cost,result.temperature);
    println!("power_coefficients_W={controls:?}");
    println!("gradient={:?} replayed_steps={} peak_checkpoints={}",result.gradient,result.replays,result.checkpoints);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn constant_power_matches_analytic_temperature_and_uniform_control_derivative() {
        let power = 20.0;
        let result = evaluate([power;8],true).unwrap();
        let response = (1.0-(-CONDUCTANCE/CAPACITY*HORIZON).exp())/CONDUCTANCE;
        let exact = AMBIENT+power*response;
        assert!((result.temperature-exact).abs()<1e-6);
        let derivative = (exact-TARGET)/(SCALE*SCALE)*response + PENALTY*2.0*power/(MAX_POWER*MAX_POWER);
        assert!((result.gradient.iter().sum::<f64>()-derivative).abs()<1e-7);
    }
    #[test]
    fn reverse_gradients_drive_bounded_control_improvement() {
        let start = evaluate([5.0;8],false).unwrap();
        let (controls,result,_) = optimize(30).unwrap();
        assert!(result.cost < 0.02*start.cost);
        assert!((result.temperature-TARGET).abs()<1.0);
        assert!(controls.iter().all(|v| *v>=0.0 && *v<=MAX_POWER));
        assert!(result.gradient.iter().all(|v|v.is_finite()));
    }
}
