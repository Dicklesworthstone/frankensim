//! Original PSE diffusion composed with the original conservative advection.

use super::*;
use crate::diffusion::moving::MovingDiffusion;
use crate::diffusion::{self, PseDiffusion};

/// Evolve viscous planar vortices with `D(dt/2) A(dt) D(dt/2)`, using the
/// original PSE diffusion and exact-pair conservative advection owners.
///
/// Positive quadrature areas move with their particles; core radius is the
/// separate common Biot–Savart regularization length. Time is nonnegative:
/// viscous diffusion is not a reversible negative-time operation. Every area
/// is admitted even for zero time, zero steps or zero viscosity. Within this
/// nonnegative-time domain, zero viscosity recovers the original conservative
/// trajectory, pair/work counts and one-vector live-memory bound exactly.
///
/// Positive time and viscosity charge `3*S*N*(N-1)` pair visits and
/// `N + S + 3*S*N*(N-1)` work units. Only one N-particle vector is owned; there
/// are no per-step allocations, RK4 buffers or additional material parameters.
/// These are operation counts, not a wall-clock performance guarantee.
///
/// Advection retains its linear/angular vortex-impulse invariants. Diffusion
/// conserves circulation and contracts its vorticity range/weighted enstrophy,
/// but does NOT exactly conserve linear/angular impulse on an arbitrary finite
/// cloud. The full viscous result inherits that limitation; it does not project
/// away physical diffusion or its finite-cloud quadrature error. Cloud/kernel
/// and temporal refinement remain distinct requirements. No finite-domain fluid
/// momentum, wall reaction, evolving core or viscous-heating claim is made.
///
/// # Errors
/// Original material/input, resource, cancellation and nonrepresentable-pair
/// refusals publish no partial cloud and leave the input untouched.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    dt: f64,
    steps: usize,
    core: f64,
    law: PseDiffusion,
    budget: VpmBudget,
) -> Result<VpmRun, VpmError> {
    diffusion::validate_inputs(particles.len(), areas_m2, dt)?;
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None));
    }
    admit_cap("particles", particles.len(), budget.max_particles)?;
    admit_cap("steps", steps, budget.max_steps)?;
    let active_advection = dt > 0.0;
    let active_diffusion = active_advection && law.kinematic_viscosity_m2_s() > 0.0;
    let sweeps = if active_diffusion { 3 } else { 1 };
    let pair_evaluations = if active_advection {
        diffusion::pair_work(particles.len(), steps, sweeps)?
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
        .and_then(|work| work.checked_add(steps as u64))
        .ok_or(VpmError::PlanOverflow {
            resource: "conservative viscous total work",
        })?;
    let peak_live_bytes = particles
        .len()
        .checked_mul(core::mem::size_of::<VortexParticle>())
        .ok_or(VpmError::PlanOverflow {
            resource: "conservative viscous live bytes",
        })?;
    admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
    let mut accountant = AdmittedBudget::admit_ambient(cx, work_units)?;
    accountant.checkpoint("fs-vpm.conservative-admission", cx)?;
    let mut state = reserved_vec("conservative-state", particles.len())?;
    let mut meter = WorkMeter::default();
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(&mut accountant, cx, "fs-vpm.conservative-input", || {
            validate_particle(particle, index)?;
            diffusion::validate_area(areas_m2[index], index)?;
            state.push(particle);
            Ok(())
        })?;
    }
    meter.boundary(&mut accountant, cx, "fs-vpm.conservative-input-complete")?;
    let mut completed_pairs = 0;
    for _ in 0..steps {
        step(
            &mut state,
            areas_m2,
            dt,
            core,
            law,
            &mut completed_pairs,
            &mut meter,
            &mut accountant,
            cx,
        )?;
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

// Shared in-place owner: scheduled remeshing must not copy the cloud or create
// a new execution account between the original D/2-A-D/2 steps.
#[allow(clippy::too_many_arguments)]
pub(super) fn step(
    state: &mut [VortexParticle],
    areas: &[f64],
    dt: f64,
    core: f64,
    law: PseDiffusion,
    completed_pairs: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    weighted_step(
        state,
        areas,
        dt,
        1.0,
        core,
        law,
        completed_pairs,
        meter,
        accountant,
        cx,
    )
}

// A separate positive fraction preserves resolvable subnormal physical time.
// The adaptive owner uses this exact existing Strang map for full/half trials.
#[allow(clippy::too_many_arguments)]
pub(super) fn weighted_step(
    state: &mut [VortexParticle],
    areas: &[f64],
    dt: f64,
    fraction: f64,
    core: f64,
    law: PseDiffusion,
    completed_pairs: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    weighted_step_with(
        state,
        areas,
        dt,
        fraction,
        core,
        law,
        &mut MovingDiffusion::dense(),
        completed_pairs,
        meter,
        accountant,
        cx,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn weighted_step_with(
    state: &mut [VortexParticle],
    areas: &[f64],
    dt: f64,
    fraction: f64,
    core: f64,
    law: PseDiffusion,
    backend: &mut MovingDiffusion,
    completed_pairs: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    let diffuse = dt > 0.0 && law.kinematic_viscosity_m2_s() > 0.0;
    if diffuse {
        backend.sweep(
            state,
            areas,
            law,
            dt,
            0.25 * fraction,
            completed_pairs,
            meter,
            accountant,
            cx,
        )?;
    }
    if dt > 0.0 {
        backend.admit_advection(state.len(), *completed_pairs, meter, accountant)?;
        super::weighted_symmetric_sweep(
            state,
            dt,
            core,
            fraction,
            completed_pairs,
            meter,
            accountant,
            cx,
        )?;
        backend.invalidate();
    }
    if diffuse {
        backend.sweep(
            state,
            areas,
            law,
            dt,
            0.25 * fraction,
            completed_pairs,
            meter,
            accountant,
            cx,
        )?;
    }
    backend.admit_work(1, meter, accountant)?;
    meter.unit(accountant, cx, "fs-vpm.conservative-step", || Ok(()))?;
    meter.boundary(accountant, cx, "fs-vpm.conservative-step-complete")
}

// Reverse positive Strang composition A(dt/2) D(dt) A(dt/2), used only as
// an independent adaptive witness. Reuse the original pair owners and keep
// half intervals as fractions so a subnormal physical dt is not lost.
#[allow(clippy::too_many_arguments)]
pub(super) fn reverse_step(
    state: &mut [VortexParticle],
    areas: &[f64],
    dt: f64,
    core: f64,
    law: PseDiffusion,
    backend: &mut MovingDiffusion,
    completed_pairs: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    if dt > 0.0 {
        backend.admit_advection(state.len(), *completed_pairs, meter, accountant)?;
        super::weighted_symmetric_sweep(
            state,
            dt,
            core,
            0.5,
            completed_pairs,
            meter,
            accountant,
            cx,
        )?;
        backend.invalidate();
        if law.kinematic_viscosity_m2_s() > 0.0 {
            backend.sweep(
                state,
                areas,
                law,
                dt,
                0.5,
                completed_pairs,
                meter,
                accountant,
                cx,
            )?;
        }
        backend.admit_advection(state.len(), *completed_pairs, meter, accountant)?;
        super::weighted_symmetric_sweep(
            state,
            dt,
            core,
            0.5,
            completed_pairs,
            meter,
            accountant,
            cx,
        )?;
        backend.invalidate();
    }
    backend.admit_work(1, meter, accountant)?;
    meter.unit(accountant, cx, "fs-vpm.viscous-witness-step", || Ok(()))?;
    meter.boundary(accountant, cx, "fs-vpm.viscous-witness-step-complete")
}
