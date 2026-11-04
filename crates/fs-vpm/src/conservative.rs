//! Momentum-preserving inviscid planar vortex evolution.
//!
//! A forward half sweep over unordered pairs followed by the exact reverse
//! half sweep is a second-order, time-reversible composition of the exact
//! [`crate::pairwise::advance_pair`] flows. Each subflow preserves both first
//! circulation moments and the quadratic circulation moment, so the composed
//! trajectory retains linear and angular impulse to floating-point roundoff.
//! No projection, circulation adjustment or global-centroid division is used.
//!
//! [`CompositionOrder::Fourth`] optionally composes three original reversible
//! sweeps for fourth-order inviscid accuracy with the same impulse invariants.
//! Finite-step phase/energy errors depend on pair order and require time
//! refinement. Conservation alone does not establish trajectory accuracy.

use super::{
    AdmittedBudget, Cx, VortexParticle, VpmBudget, VpmError, VpmRun, WorkMeter, admit_cap,
    invalid_scalar, pairwise::advance_validated_pair, reserved_vec, validate_particle,
};

mod viscous;
pub use viscous::simulate_with_cx as simulate_viscous_with_cx;

/// Adaptive positive-time inviscid/source and moving viscous particle windows.
pub mod adaptive;

mod composition;
use composition::Composition;
pub use composition::CompositionOrder;

/// Checked viscous evolution over an authored sequence of regular remeshes.
pub mod remeshed;

/// Prescribed incompressible affine flow around the original native evolution.
pub mod ambient;

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
    simulate_with_order_cx(
        cx,
        particles,
        dt,
        steps,
        core,
        CompositionOrder::Second,
        budget,
    )
}

/// Choose second- or fourth-order inviscid evolution in the original transaction.
///
/// Fourth order uses `S(w1*dt) S(w0*dt) S(w1*dt)` with the original symmetric
/// sweep `S`, `w1=1/(2-cbrt(2))` and `w0=1-2*w1`. Each nonzero physical step
/// costs `3*N*(N-1)` pair visits and one step unit. The same single particle
/// vector is reused; there are no stage or per-step allocations. Negative
/// internal intervals are reversible vortex evolution, not negative viscosity.
///
/// The original entry point and [`CompositionOrder::Second`] are identical,
/// including work counts and checkpoint order. Zero-time/step and empty or
/// single-particle runs retain their original work counts at either order.
/// Smooth finite-core/noncolliding evolution admits the fourth-order time
/// refinement claim; singular collisions, exact many-vortex energy, finite-volume
/// momentum, adaptive error control and wall-clock speed are not claimed.
///
/// # Errors
/// Original input/resource/numerical/cancellation refusals, including any
/// unrepresentable internal pair update. No partial state escapes.
pub fn simulate_with_order_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    dt: f64,
    steps: usize,
    core: f64,
    order: CompositionOrder,
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
        super::diffusion::pair_work(particles.len(), steps, order.sweeps() as u64)?
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
    let composition = Composition::new(order);
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
            for &weight in composition.weights() {
                weighted_symmetric_sweep(
                    &mut state,
                    dt,
                    core,
                    weight,
                    &mut completed_pairs,
                    &mut meter,
                    &mut accountant,
                    cx,
                )?;
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

// The inviscid and viscous transactions share exactly the same pair map,
// ordering, arithmetic and work-meter phases. This performs one full advection
// interval; diffusion composes around it without replacing the advection owner.
fn symmetric_sweep(
    state: &mut [VortexParticle],
    dt: f64,
    core: f64,
    completed_pairs: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    weighted_symmetric_sweep(state, dt, core, 1.0, completed_pairs, meter, accountant, cx)
}

// The existing pair leaf carries a positive time fraction in its scaled
// denominator. Put the sign in dt, preserving subnormal and extreme products
// without exposing negative fractions to that leaf's magnitude arithmetic.
#[allow(clippy::too_many_arguments)]
fn weighted_symmetric_sweep(
    state: &mut [VortexParticle],
    dt: f64,
    core: f64,
    weight: f64,
    completed_pairs: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    let signed_dt = if weight.is_sign_negative() { -dt } else { dt };
    let half_fraction = 0.5 * weight.abs();
    for i in 0..state.len() {
        for j in (i + 1)..state.len() {
            meter.unit(accountant, cx, "fs-vpm.conservative-forward", || {
                advance(state, i, j, signed_dt, core, half_fraction)?;
                *completed_pairs += 1;
                Ok(())
            })?;
        }
    }
    for i in (0..state.len()).rev() {
        for j in ((i + 1)..state.len()).rev() {
            meter.unit(accountant, cx, "fs-vpm.conservative-reverse", || {
                advance(state, i, j, signed_dt, core, half_fraction)?;
                *completed_pairs += 1;
                Ok(())
            })?;
        }
    }
    Ok(())
}

fn advance(
    state: &mut [VortexParticle],
    i: usize,
    j: usize,
    dt: f64,
    core: f64,
    half_fraction: f64,
) -> Result<(), VpmError> {
    let pair =
        advance_validated_pair([state[i], state[j]], dt, core, half_fraction).map_err(|_| {
            VpmError::NonFiniteResult {
                phase: "fs-vpm.conservative-pair",
                particle: i,
            }
        })?;
    state[i] = pair[0];
    state[j] = pair[1];
    Ok(())
}
