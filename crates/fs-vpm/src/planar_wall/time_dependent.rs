//! Chronological material-viscosity exposure at a stationary stress-free wall.

use super::PlanarWall2;
use crate::conservative::{TimeDependentVpmError, TimeDependentVpmRun};
use crate::diffusion;
use crate::diffusion::viscosity_history::PseViscosityHistory;
use crate::{
    AdmittedBudget, RkScratch, VelocityKernel, VortexParticle, VpmBudget, VpmError, VpmPlan,
    VpmRun, WorkMeter, admit_cap, checked_rk4_step, invalid_scalar, reserved_vec,
    validate_particle,
};
use fs_exec::{BudgetRefusal, Cx};

/// Advect and diffuse at a flat stress-free wall under a homogeneous `nu(t)`.
///
/// The original viscosity history defines the coefficient law, including its
/// explicit interpolation and source support. Each interval ends at the next
/// coefficient knot, the requested endpoint or the maximum represented step.
/// The first and second physical halves use their own exactly integrated mean
/// viscosities around the original real-plus-image RK4 advection. A full-step
/// mean on both sides would misorder wall loss and the motion induced by the
/// changing particle strengths.
///
/// Odd-image PSE imposes zero bulk vorticity at the stationary flat wall.
/// Signed bulk circulation can leave through that boundary; images are never
/// physical particles or an additional circulation inventory. This is neither
/// no-slip vorticity generation nor a stored boundary vortex sheet or wall load.
/// Spatially varying viscosity, changing density and thermal feedback are not
/// inferred from a prescribed homogeneous coefficient history.
///
/// The composition is second order on smooth history segments and requires
/// independent time and cloud/kernel refinement. Constant histories on the
/// same represented step partition retain the original constant-law particle
/// bits and pair/work/memory receipt. Each step charges `8*N*N` RK4 pair
/// visits, plus `2*N*N` for each positive-viscosity half. Work additionally
/// includes `N` input visits and `4*N+1` visits per step.
///
/// One original RK4 scratch set is reused, with no per-step allocation, timeline
/// vector or material callback. Planning is bounded by the step cap and polls
/// cancellation every 256 intervals before one aggregate ambient admission.
/// No particle state is returned until the whole physical window succeeds.
/// Retain any `MaterialPseViscosityHistory` binding for its original source
/// owners and transport snapshots, just as for its free-space consumer.
///
/// # Errors
/// Refuses unsupported coefficient windows, invalid or stagnant time steps,
/// original numerical inputs, any particle/stage crossing the wall, resource
/// limits, cancellation, allocation failure or nonrepresentable arithmetic.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    start_s: f64,
    end_s: f64,
    maximum_step_s: f64,
    core: f64,
    history: PseViscosityHistory<'_>,
    wall: PlanarWall2,
    budget: VpmBudget,
) -> Result<TimeDependentVpmRun, TimeDependentVpmError> {
    diffusion::validate_inputs(particles.len(), areas_m2, end_s - start_s)?;
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None).into());
    }
    admit_cap("particles", particles.len(), budget.max_particles)?;
    let profile = history.profile();
    let mut time = start_s;
    let mut steps = 0_usize;
    let mut diffusion_visits = 0_u64;
    let mut pair_evaluations = 0_u64;
    while let Some(interval) = profile.next_interval(time, end_s, maximum_step_s)? {
        if steps % 256 == 0 {
            cx.checkpoint().map_err(|_| BudgetRefusal::Cancelled {
                phase: "fs-vpm.wall-viscosity-history-plan",
            })?;
        }
        steps = steps.checked_add(1).ok_or(VpmError::PlanOverflow {
            resource: "wall viscosity history steps",
        })?;
        admit_cap("steps", steps, budget.max_steps)?;
        let positive_halves = interval
            .half_mean_viscosities_m2_s
            .iter()
            .filter(|&&nu| nu > 0.0)
            .count() as u64;
        diffusion_visits = diffusion_visits
            .checked_add(diffusion::boundary_pair_work(
                particles.len(),
                1,
                positive_halves,
                Some(wall),
            )?)
            .ok_or(VpmError::PlanOverflow {
                resource: "wall viscosity history diffusion work",
            })?;
        pair_evaluations = VpmPlan::checked(particles.len(), steps)?
            .pair_evaluations
            .checked_mul(2)
            .and_then(|visits| visits.checked_add(diffusion_visits))
            .ok_or(VpmError::PlanOverflow {
                resource: "wall viscosity history pair work",
            })?;
        if pair_evaluations > budget.max_pair_evaluations {
            return Err(VpmError::CapExceeded {
                resource: "pair evaluations",
                required: pair_evaluations.into(),
                maximum: budget.max_pair_evaluations.into(),
            }
            .into());
        }
        time = interval.time_s[1];
    }
    let mut plan = VpmPlan::checked(particles.len(), steps)?;
    // Images add one original RK4 pair plan. Diffusion has already counted
    // only active halves, including physical self-image decay at N=1.
    let additional_visits =
        plan.pair_evaluations
            .checked_add(diffusion_visits)
            .ok_or(VpmError::PlanOverflow {
                resource: "wall viscosity history additional work",
            })?;
    plan.work_units =
        plan.work_units
            .checked_add(additional_visits)
            .ok_or(VpmError::PlanOverflow {
                resource: "wall viscosity history total work",
            })?;
    plan.pair_evaluations = pair_evaluations;
    admit_cap(
        "logical live bytes",
        plan.peak_live_bytes,
        budget.max_live_bytes,
    )?;
    let mut accountant = AdmittedBudget::admit_ambient(cx, plan.work_units)?;
    accountant.checkpoint("fs-vpm.wall-viscosity-history-admission", cx)?;
    let mut state = reserved_vec("wall-viscosity-history-state", particles.len())?;
    let mut scratch = RkScratch::try_new(if steps == 0 { 0 } else { particles.len() })?;
    let mut meter = WorkMeter::default();
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(
            &mut accountant,
            cx,
            "fs-vpm.wall-viscosity-history-input",
            || {
                validate_particle(particle, index)?;
                diffusion::validate_area(areas_m2[index], index)?;
                wall.admit_position(particle.pos, index)?;
                state.push(particle);
                Ok(())
            },
        )?;
    }
    meter.boundary(
        &mut accountant,
        cx,
        "fs-vpm.wall-viscosity-history-input-complete",
    )?;
    let kernel = VelocityKernel::new(core);
    let mut completed_pairs = 0;
    time = start_s;
    while let Some(interval) = profile.next_interval(time, end_s, maximum_step_s)? {
        let [left, right] = interval.time_s;
        let dt = right - left;
        let [first, second] = history.half_laws(interval)?;
        if first.kinematic_viscosity_m2_s() > 0.0 {
            diffusion::boundary_sweep(
                &mut state,
                areas_m2,
                first,
                dt,
                0.25,
                &mut completed_pairs,
                &mut meter,
                &mut accountant,
                cx,
                Some(wall),
            )?;
        }
        checked_rk4_step(
            &mut state,
            dt,
            kernel,
            Some(wall),
            &mut scratch,
            &mut completed_pairs,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        if second.kinematic_viscosity_m2_s() > 0.0 {
            diffusion::boundary_sweep(
                &mut state,
                areas_m2,
                second,
                dt,
                0.25,
                &mut completed_pairs,
                &mut meter,
                &mut accountant,
                cx,
                Some(wall),
            )?;
        }
        meter.boundary(&mut accountant, cx, "fs-vpm.step-complete")?;
        time = right;
    }
    meter.flush(&mut accountant, "fs-vpm.wall-viscosity-history-finalize")?;
    let consumption = accountant.consumption();
    if meter.completed != plan.work_units
        || completed_pairs != plan.pair_evaluations
        || consumption.cost_charged != plan.work_units
    {
        return Err(VpmError::AccountingMismatch {
            planned: plan.work_units,
            completed: meter.completed,
            charged: consumption.cost_charged,
        }
        .into());
    }
    accountant.checkpoint("fs-vpm.wall-viscosity-history-publication", cx)?;
    Ok(TimeDependentVpmRun {
        run: VpmRun {
            particles: state,
            steps_completed: steps,
            pair_evaluations: completed_pairs,
            work_units: plan.work_units,
            peak_live_bytes: plan.peak_live_bytes,
            budget: accountant.consumption(),
        },
        time_s: [start_s, time],
    })
}
