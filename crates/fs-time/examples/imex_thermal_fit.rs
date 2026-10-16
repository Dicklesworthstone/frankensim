//! Fit conductance and heater power through an implicit thermal trajectory.
//!
//! `cargo run -p fs-time --example imex_thermal_fit`
//!
//! Two lumped bodies exchange heat and cool to a zero-reference reservoir.
//! Temperatures are excess kelvin, capacities are 1 and 2 J/K, conductances
//! are W/K, heater power is W, and time is seconds. Synthetic endpoint data
//! come from this same discrete model: this demonstrates inverse computation,
//! not independent physical validation or parameter identifiability in general.
//! Forty-step gradients use at most six parked checkpoints with bounded replay.
use fs_solver::LinearOp;
use fs_time::stiff::{
    IdentityPreconditioner, ImexSolveConfig, OperatorImex2,
    adjoint::ImexVjp,
    adjoint::trajectory::{
        ImexRecordingConfig, ImexRecordingStatus, ImexReplayBudget, RecordedImex2,
    },
};

struct Thermal {
    conductance: f64,
    heater: f64,
}
impl LinearOp for Thermal {
    fn n(&self) -> usize {
        2
    }
    fn apply(&self, u: &[f64], out: &mut [f64]) {
        let exchange = self.conductance * (u[1] - u[0]);
        out[0] = exchange - 0.4 * u[0];
        out[1] = 0.5 * (-exchange - 0.7 * u[1]);
    }
    fn apply_transpose(&self, u: &[f64], out: &mut [f64]) {
        out[0] = -(self.conductance + 0.4) * u[0] + 0.5 * self.conductance * u[1];
        out[1] = self.conductance * u[0] - 0.5 * (self.conductance + 0.7) * u[1];
    }
}
impl ImexVjp for Thermal {
    fn parameter_count(&self) -> usize {
        2
    }
    fn nonlinear(&self, _u: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[self.heater, 0.0]);
    }
    fn nonlinear_vjp(
        &self,
        _u: &[f64],
        seed: &[f64],
        ub: &mut [f64],
        pb: &mut [f64],
    ) -> Result<(), String> {
        ub.fill(0.0);
        pb.copy_from_slice(&[0.0, seed[0]]);
        Ok(())
    }
    fn linear_parameter_vjp(&self, u: &[f64], seed: &[f64], pb: &mut [f64]) -> Result<(), String> {
        pb.copy_from_slice(&[(u[1] - u[0]) * (seed[0] - 0.5 * seed[1]), 0.0]);
        Ok(())
    }
}

fn forward<'a>(
    method: &OperatorImex2,
    model: &'a Thermal,
) -> Result<RecordedImex2<'a, Thermal, IdentityPreconditioner>, Box<dyn std::error::Error>> {
    let config = ImexRecordingConfig {
        steps: 40,
        max_workspace_components: method
            .adjoint_workspace_components(2)
            .ok_or("workspace overflow")?,
    };
    let mut tape = RecordedImex2::new(
        *method,
        model,
        &IdentityPreconditioner,
        0.0,
        &[1.0, 0.2],
        config,
    )?;
    if tape.advance(40, 40, &mut || false)?.status != ImexRecordingStatus::ReachedEnd {
        return Err("forward trajectory incomplete".into());
    }
    Ok(tape)
}

fn loss(endpoint: &[f64], target: &[f64]) -> f64 {
    0.5 * endpoint
        .iter()
        .zip(target)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let method = OperatorImex2::new(
        2,
        0.05,
        ImexSolveConfig {
            tolerance: 1e-13,
            restart: 2,
            max_cycles: 4,
        },
    );
    let target = forward(
        &method,
        &Thermal {
            conductance: 1.2,
            heater: 0.8,
        },
    )?
    .state()
    .to_vec();
    let mut model = Thermal {
        conductance: 0.5,
        heater: 0.3,
    };
    println!("iteration,loss_k2,conductance_w_per_k,heater_w");
    for iteration in 0..=80 {
        let tape = forward(&method, &model)?;
        let value = loss(tape.state(), &target);
        if iteration % 10 == 0 {
            println!(
                "{iteration},{value:.12e},{:.9},{:.9}",
                model.conductance, model.heater
            );
        }
        if iteration == 80 || value < 1e-14 {
            break;
        }
        let bar = tape
            .state()
            .iter()
            .zip(&target)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        let gradient = tape
            .pullback(
                &bar,
                &[0.0; 2],
                &IdentityPreconditioner,
                ImexReplayBudget {
                    checkpoints: 6,
                    forward_steps: 200,
                },
                &mut || false,
            )?
            .parameters;
        let mut step = 8.0;
        let mut accepted = false;
        for _ in 0..20 {
            let candidate = Thermal {
                conductance: model.conductance - step * gradient[0],
                heater: model.heater - step * gradient[1],
            };
            if candidate.conductance > 0.0 && candidate.heater >= 0.0 {
                if loss(forward(&method, &candidate)?.state(), &target)
                    < value - 1e-4 * step * gradient.iter().map(|v| v * v).sum::<f64>()
                {
                    model = candidate;
                    accepted = true;
                    break;
                }
            }
            step *= 0.5;
        }
        if !accepted {
            return Err("bounded line search exhausted before convergence".into());
        }
    }
    Ok(())
}
