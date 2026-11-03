//! Momentum-preserving inviscid planar vortex evolution.
//!
//! A forward half sweep over unordered pairs followed by the exact reverse
//! half sweep is a second-order, time-reversible composition of the exact
//! [`crate::pairwise::advance_pair`] flows. Each subflow preserves both first
//! circulation moments and the quadratic circulation moment, so the composed
//! trajectory retains linear and angular impulse to floating-point roundoff.
//! No projection, circulation adjustment or global-centroid division is used.
//!
//! This is a separate choice from the existing fourth-order RK4 entry. Its
//! finite-step phase/energy errors depend on pair order and require time
//! refinement. Conservation alone does not establish trajectory accuracy.

use super::{
    AdmittedBudget, Cx, VortexParticle, VpmBudget, VpmError, VpmRun, WorkMeter, admit_cap,
    invalid_scalar, pairwise::advance_validated_pair, reserved_vec, validate_particle,
};

/// Advance the inviscid cloud by symmetric exact-pair splitting under the
/// existing local and ambient execution budgets.
///
/// Uses the same positive-counterclockwise algebraic-core Biot–Savart model
/// as [`crate::simulate_with_cx`]. Negative time is supported. Each nonzero
/// step visits `N(N-1)` pairs, half in lexicographic order and half in its exact
/// reverse. Zero time skips pair work. Every step also costs one boundary unit;
/// input validation/copy costs `N`. Only one N-particle payload is allocated,
/// with no per-step or pair allocations. No partial state escapes on refusal.
///
/// Exact pair rotation preserves `sum Gamma*x`, `sum Gamma*y`, and
/// `sum Gamma*|x|²`, including mixed signs and zero total circulation. These
/// point-vortex impulse moments do not assert finite-volume fluid momentum,
/// wall reaction forces, viscosity, or exact many-vortex energy conservation.
///
/// # Errors
/// Refuses invalid inputs, overflow of the checked work plan, local or ambient
/// budget exhaustion, cancellation/deadline, allocation failure, and any pair
/// whose separation, rotation angle or final position is not representable.
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    budget: VpmBudget,
) -> Result<VpmRun, VpmError> {
    if !dt.is_finite() {
        return Err(invalid_scalar("dt", dt, None));
    }
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None));
    }
    admit_cap("particles", particles.len(), budget.max_particles)?;
    admit_cap("steps", steps, budget.max_steps)?;
    let active = dt != 0.0;
    let pair_evaluations = if active {
        super::diffusion::pair_work(particles.len(), steps, 1)?
    } else {
        0
    };
    if pair_evaluations > budget.max_pair_evaluations {
        return Err(VpmError::CapExceeded {
            resource: "pair evaluations",
            required: pair_evaluations.into(),
            maximum: budget.max_pair_evaluations.into(),
        });
    }
    let work_units = pair_evaluations
        .checked_add(particles.len() as u64)
        .and_then(|value| value.checked_add(steps as u64))
        .ok_or(VpmError::PlanOverflow {
            resource: "conservative total work",
        })?;
    let peak_live_bytes = particles
        .len()
        .checked_mul(core::mem::size_of::<VortexParticle>())
        .ok_or(VpmError::PlanOverflow {
            resource: "conservative live bytes",
        })?;
    admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
    let mut accountant = AdmittedBudget::admit_ambient(cx, work_units)?;
    accountant.checkpoint("fs-vpm.conservative-admission", cx)?;
    let mut state = reserved_vec("conservative-state", particles.len())?;
    let mut meter = WorkMeter::default();
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(&mut accountant, cx, "fs-vpm.conservative-input", || {
            validate_particle(particle, index)?;
            state.push(particle);
            Ok(())
        })?;
    }
    meter.boundary(&mut accountant, cx, "fs-vpm.conservative-input-complete")?;
    let mut completed_pairs = 0;
    for _ in 0..steps {
        if active {
            for i in 0..state.len() {
                for j in (i + 1)..state.len() {
                    meter.unit(&mut accountant, cx, "fs-vpm.conservative-forward", || {
                        advance(&mut state, i, j, dt, core)?;
                        completed_pairs += 1;
                        Ok(())
                    })?;
                }
            }
            for i in (0..state.len()).rev() {
                for j in ((i + 1)..state.len()).rev() {
                    meter.unit(&mut accountant, cx, "fs-vpm.conservative-reverse", || {
                        advance(&mut state, i, j, dt, core)?;
                        completed_pairs += 1;
                        Ok(())
                    })?;
                }
            }
        }
        meter.unit(&mut accountant, cx, "fs-vpm.conservative-step", || Ok(()))?;
        meter.boundary(&mut accountant, cx, "fs-vpm.conservative-step-complete")?;
    }
    meter.flush(&mut accountant, "fs-vpm.conservative-finalize")?;
    let consumption = accountant.consumption();
    if meter.completed != work_units
        || completed_pairs != pair_evaluations
        || consumption.cost_charged != work_units
    {
        return Err(VpmError::AccountingMismatch {
            planned: work_units,
            completed: meter.completed,
            charged: consumption.cost_charged,
        });
    }
    accountant.checkpoint("fs-vpm.conservative-publication", cx)?;
    Ok(VpmRun {
        particles: state,
        steps_completed: steps,
        pair_evaluations,
        work_units,
        peak_live_bytes,
        budget: accountant.consumption(),
    })
}

fn advance(
    state: &mut [VortexParticle],
    i: usize,
    j: usize,
    dt: f64,
    core: f64,
) -> Result<(), VpmError> {
    let pair = advance_validated_pair([state[i], state[j]], dt, core, 0.5).map_err(|_| {
        VpmError::NonFiniteResult {
            phase: "fs-vpm.conservative-pair",
            particle: i,
        }
    })?;
    state[i] = pair[0];
    state[j] = pair[1];
    Ok(())
}
