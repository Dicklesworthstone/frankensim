//! Static loading of the original conservative constitutive potentials.
//!
//! One factorization of the original bilateral stiffness preconditions descent
//! of U_linear(q) + V(q) - F.q. A positive secant curvature predicts the next
//! step length; backtracking uses the same stored energy and
//! exact gradient as playback, including mapped replacement stiffness. There
//! is no settling time integration, invented damping, retuning or modal reset.
//! A stationary result is accepted only after every actual generalized force
//! closes its declared allowance. This is not a uniqueness or stability proof.
use super::*;
use crate::render::schedule::force::coupled::equilibrium::StaticResponse;

/// A force-loaded initial condition, before any physical sample is rendered.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearPreloadDiagnostics {
    /// Constitutive trial evaluations, including rejected backtracking trials.
    /// Bounded by the original `NonlinearModalConfig::max_iterations`.
    pub iterations: usize,
    /// Absolute force residual at the largest residual/allowance ratio
    /// [N/sqrt(kg)], not a displacement convergence criterion.
    pub force_residual: f64,
    /// Allowance at that same coordinate [N/sqrt(kg)]. Every coordinate passed.
    pub force_tolerance: f64,
    /// Conservative correction to the original linear storage [J]. May be
    /// signed when the full potential replaces selected diagonal stiffness.
    pub nonlinear_energy_j: f64,
    /// Complete stored energy [J]. This is not pre-window actuator work.
    pub total_energy_j: f64,
}

impl FrictionModalSystem {
    /// Last accepted nonlinear static initialization, without advancing time.
    pub fn nonlinear_preload_diagnostics(&self) -> Option<&NonlinearPreloadDiagnostics> {
        self.nonlinear
            .as_ref()
            .and_then(|state| state.preload.as_ref())
    }
}

struct Evaluation {
    energy_j: f64,
    objective_j: f64,
    work_scale_j: f64,
    gradient: Vec<f64>,
    reaction: Vec<f64>,
}

// The static diagonal and bilateral springs are exactly the original network
// stiffness. All hereditary branches are relaxed; dashpots carry no static
// force. Constitutive energy/gradient remain with their original Storage owner.
fn evaluate(
    state: &NonlinearFrictionState,
    network: &CoupledModalSystem,
    q: &[f64],
    external: &[f64],
    gate: &CancelGate,
) -> Result<Evaluation, ModalFrictionError> {
    poll(Some(gate))?;
    let mut energy = state.potential_energy(q)?;
    let mut reaction = vec![0.0; q.len()];
    state.reactions(q, q, &mut reaction, Some(gate))?;
    let mut gradient = external
        .iter()
        .zip(&reaction)
        .map(|(force, reaction)| finite(-force - reaction))
        .collect::<Result<Vec<_>, _>>()?;
    let mut index = 0;
    for model in &network.models {
        poll(Some(gate))?;
        for mode in model.modes() {
            let force =
                finite(mode.angular_frequency_rad_s * mode.angular_frequency_rad_s * q[index])?;
            gradient[index] = finite(gradient[index] + force)?;
            energy = finite(energy + 0.5 * force * q[index])?;
            index += 1;
        }
    }
    for (column, link) in network.columns.iter().zip(&network.connections) {
        poll(Some(gate))?;
        let extension = finite(dot(column, q)? - link.rest_extension_m)?;
        let force = finite(link.stiffness_n_m * extension)?;
        energy = finite(energy + 0.5 * force * extension)?;
        for (gradient, shape) in gradient.iter_mut().zip(column) {
            *gradient = finite(*gradient + shape * force)?;
        }
    }
    let work = dot(external, q)?;
    let work_scale = external
        .iter()
        .zip(q)
        .try_fold(0.0, |sum, (f, q)| finite(sum + (f * q).abs()))?;
    Ok(Evaluation {
        energy_j: energy,
        objective_j: finite(energy - work)?,
        work_scale_j: work_scale,
        gradient,
        reaction,
    })
}

struct Residual {
    allowance: Vec<f64>,
    force: f64,
    tolerance: f64,
    closed: bool,
}

// A spectral step in the original stiffness metric. For s=q_new-q_old,
// alpha=(s^T K s)/(s^T (grad_new-grad_old)). Compute its numerator directly
// from positive spring squares, avoiding cancellation of large nonlinear
// forces. A nonpositive/unrepresentable secant falls back to unit step;
// every prediction still passes the same energy and force safeguards.
fn next_step(
    network: &CoupledModalSystem,
    before: &[f64],
    after: &[f64],
    old: &Evaluation,
    new: &Evaluation,
    gate: &CancelGate,
) -> Result<f64, ModalCouplingError> {
    let s: Vec<_> = before.iter().zip(after).map(|(a, b)| b - a).collect();
    let y: Vec<_> = old
        .gradient
        .iter()
        .zip(&new.gradient)
        .map(|(a, b)| b - a)
        .collect();
    let denominator = dot(&s, &y)?;
    let mut numerator = 0.0;
    for (mode, displacement) in network
        .models
        .iter()
        .flat_map(|model| model.modes())
        .zip(&s)
    {
        let scaled = finite(mode.angular_frequency_rad_s * displacement)?;
        numerator = finite(numerator + scaled * scaled)?;
    }
    for (column, link) in network.columns.iter().zip(&network.connections) {
        poll(Some(gate))?;
        let extension = dot(column, &s)?;
        numerator = finite(numerator + link.stiffness_n_m * extension * extension)?;
    }
    let predicted = numerator / denominator;
    Ok(
        if numerator > 0.0 && denominator > 0.0 && predicted.is_finite() && predicted > 0.0 {
            predicted
        } else {
            1.0
        },
    )
}
fn residual(
    value: &Evaluation,
    external: &[f64],
    config: NonlinearModalConfig,
) -> Result<Residual, ModalCouplingError> {
    let mut allowance = Vec::with_capacity(external.len());
    let (mut worst, mut force, mut tolerance) = (0.0_f64, 0.0, config.force_absolute_tolerance);
    for ((gradient, reaction), external) in value.gradient.iter().zip(&value.reaction).zip(external)
    {
        let linear = finite(gradient + external + reaction)?;
        let scale = external.abs().max(reaction.abs()).max(linear.abs());
        let allowed =
            finite(config.force_absolute_tolerance + config.force_relative_tolerance * scale)?;
        allowance.push(allowed);
        let ratio = finite(gradient.abs() / allowed)?;
        if ratio > worst {
            worst = ratio;
            force = gradient.abs();
            tolerance = allowed;
        }
    }
    Ok(Residual {
        allowance,
        force,
        tolerance,
        closed: worst <= 1.0,
    })
}

pub(super) fn initialize(
    state: &NonlinearFrictionState,
    network: &mut CoupledModalSystem,
    external: &[f64],
    gate: &CancelGate,
) -> Result<NonlinearPreloadDiagnostics, ModalFrictionError> {
    poll(Some(gate))?;
    network.require_static_initialization(external)?;
    let response = StaticResponse::new(network, gate)?;
    let mut q = vec![0.0; network.mode_count()];
    let mut value = evaluate(state, network, &q, external, gate)?;
    let mut iterations = 0;
    let mut predicted_step = 1.0;
    loop {
        poll(Some(gate))?;
        let check = residual(&value, external, state.config)?;
        if check.closed {
            let forces = external
                .iter()
                .zip(&value.reaction)
                .map(|(f, r)| finite(f + r))
                .collect::<Result<Vec<_>, _>>()?;
            // Final physical caps and relaxed Maxwell history use the original
            // transactional static owner. No state is published by a trial.
            let linear = network.stage_static_equilibrium_with_allowance(
                &q,
                &forces,
                Some(&check.allowance),
                gate,
            )?;
            let nonlinear = state.potential_energy(&q)?;
            let energy = finite(linear + nonlinear)?;
            if energy < 0.0 {
                return Err(invalid("nonlinear static total storage is negative").into());
            }
            limit(
                "nonlinear static total energy",
                energy,
                network.config.maximum_total_energy_j,
            )?;
            poll(Some(gate))?;
            std::mem::swap(&mut network.models, &mut network.candidates);
            std::mem::swap(
                &mut network.material_memory,
                &mut network.candidate_material_memory,
            );
            return Ok(NonlinearPreloadDiagnostics {
                iterations,
                force_residual: check.force,
                force_tolerance: check.tolerance,
                nonlinear_energy_j: nonlinear,
                total_energy_j: energy,
            });
        }
        let failed = |iterations| ModalFrictionError::NonlinearEquilibrium {
            residual: check.force,
            tolerance: check.tolerance,
            iterations,
        };
        if iterations == state.config.max_iterations {
            return Err(failed(iterations));
        }
        let negative_gradient: Vec<_> = value.gradient.iter().map(|g| -g).collect();
        // The static response admits exactly its existing supported free
        // coordinates and unloaded gauges. Nonlinear-only support is not added.
        let direction = response.solve(network, &negative_gradient, false, false, gate)?;
        let slope = dot(&value.gradient, &direction)?;
        if slope >= 0.0 {
            return Err(failed(iterations));
        }
        let mut step = predicted_step;
        loop {
            if iterations == state.config.max_iterations {
                return Err(failed(iterations));
            }
            poll(Some(gate))?;
            let trial_q: Vec<_> = q
                .iter()
                .zip(&direction)
                .map(|(q, d)| q + step * d)
                .collect();
            if trial_q == q {
                return Err(failed(iterations));
            }
            iterations += 1;
            match evaluate(state, network, &trial_q, external, gate) {
                Ok(trial) => {
                    // Armijo descent, with only arithmetic roundoff admitted.
                    // This allowance never replaces the physical force gate.
                    let roundoff = 32.0
                        * f64::EPSILON
                        * (value.energy_j.abs()
                            + trial.energy_j.abs()
                            + value.work_scale_j
                            + trial.work_scale_j);
                    if trial.objective_j <= value.objective_j + 1e-4 * step * slope + roundoff {
                        predicted_step = next_step(network, &q, &trial_q, &value, &trial, gate)?;
                        q = trial_q;
                        value = trial;
                        break;
                    }
                }
                Err(ModalFrictionError::Network(ModalCouplingError::Cancelled)) => {
                    return Err(ModalCouplingError::Cancelled.into());
                }
                // A large trial may leave a source strain/force domain. Shrink
                // it; a valid stationary endpoint is still required to publish.
                Err(_) => {}
            }
            step *= 0.5;
        }
    }
}
