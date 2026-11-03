//! Original PSE diffusion composed with the original conservative advection.

use super::*;
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
        if active_diffusion {
            diffusion::symmetric_sweep(
                &mut state,
                areas_m2,
                law,
                dt,
                0.25,
                &mut completed_pairs,
                &mut meter,
                &mut accountant,
                cx,
            )?;
        }
        if active_advection {
            super::symmetric_sweep(
                &mut state,
                dt,
                core,
                &mut completed_pairs,
                &mut meter,
                &mut accountant,
                cx,
            )?;
        }
        if active_diffusion {
            diffusion::symmetric_sweep(
                &mut state,
                areas_m2,
                law,
                dt,
                0.25,
                &mut completed_pairs,
                &mut meter,
                &mut accountant,
                cx,
            )?;
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
