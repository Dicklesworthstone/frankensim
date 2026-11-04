//! Adaptive positive-time advection and diffusion of planar vortices.
//!
//! [`simulate_with_cx`] selects second/fourth-order inviscid pair maps, and
//! [`simulate_with_affine_flow_cx`] includes one held affine source. These
//! impulse-preserving paths share the original adaptive execution vocabulary.
//!
//! The viscous [`simulate_viscous_with_cx`] path retains positive diffusion:
//! One original D(h/2) A(h) D(h/2) step is compared with two half steps.
//! Positive viscosity also checks an independent A(h/2) D(h) A(h/2) witness,
//! exposing initial motion that fast diffusion can erase in both trial states.
//! Positions and circulation must satisfy their own physical tolerance;
//! accepting the fine state preserves the original positive diffusion maps.
//! Areas remain attached to their particles. Each diffusion stage evaluates
//! the current geometry, so no fixed-site neighbor graph is reused after motion.
//!
//! This controls temporal splitting, not spatial quadrature or fluid-model
//! error. The original direct O(N^2) pair owners remain authoritative. A
//! changing cloud still needs sufficient kernel overlap and explicit remeshing.

use crate::diffusion::adaptive::{
    AdaptiveOptions as DiffusionAdaptiveOptions, admit_count, count, next_step, overflow,
    scaled_error, trial_time, validate_options,
};
use crate::diffusion::{self, PseDiffusion};
use crate::{
    AdmittedBudget, Cx, VortexParticle, VpmError, WorkMeter, admit_cap, invalid_scalar,
    reserved_vec, validate_particle,
};
use fs_math::det;

pub use crate::diffusion::adaptive::{AdaptiveBudget, AdaptiveError, AdaptiveRun};

mod inviscid;
pub use inviscid::{
    AdaptiveOptions, AmbientAdaptiveBudget, AmbientAdaptiveRun, simulate_with_affine_flow_cx,
    simulate_with_cx,
};

/// Independent physical tolerances for vorticity and particle trajectories.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousAdaptiveOptions {
    /// Positive step limits, absolute vorticity [1/s] and relative circulation
    /// tolerances, shared with fixed-site adaptive diffusion.
    pub time: DiffusionAdaptiveOptions,
    /// Strictly positive absolute Euclidean position tolerance [m].
    /// No tolerance is scaled by distance from an arbitrary coordinate origin.
    pub absolute_position_tolerance_m: f64,
}

/// Evolve moving viscous particles for one complete nonnegative duration.
///
/// Every trial compares the same original second-order viscous composition
/// at full and half steps. The maximum of the circulation/vorticity estimate
/// and ||x_fine-x_coarse||/(3*position_tolerance) must be at most one.
/// Positive viscosity additionally compares the fine state with an independent
/// A(h/2) D(h) A(h/2) witness, without the step-doubling division by three.
/// This exposes initial transport that rapid diffusion can suppress in both
/// ordinary candidates. All estimates must pass; they remain local estimators,
/// not a certificate for arbitrary stiff transients or global trajectory error.
/// Only the two-half-step state is accepted, with no Richardson projection,
/// negative diffusion interval or post-step circulation/impulse correction.
///
/// With P=N(N-1), an active viscous trial uses 12P pair visits, of which
/// 7P evaluate diffusion kernels at the current substage positions. Its work
/// is 4N + pair_visits + 5 units: initial trial copies/estimates, witness
/// copies/estimates, four step boundaries and one trial boundary. Zero
/// viscosity skips the witness, using 3P advection visits, no diffusion
/// evaluations and 2N + pair_visits + 4 units. Initial validation costs N.
/// Three particle buffers are reserved once; the witness reuses the coarse
/// buffer. Runs with zero duration or fewer than two particles preserve the
/// inventory exactly, using one N-particle buffer and N validation units.
/// All areas and options are still checked.
///
/// The original exact pair advection preserves vortex impulse. PSE conserves
/// circulation and contracts sampled-vorticity range/enstrophy but does not
/// exactly preserve impulse on an arbitrary finite cloud. This adaptive owner
/// retains those distinct properties and does not infer wall reactions, viscous
/// heating, a changing velocity core, spatial adaptivity or a global error bound.
///
/// # Errors
/// Original physical/numerical refusal, invalid tolerances, cancellation,
/// aggregate pair/work/attempt/storage caps, an unreachable tolerance floor,
/// or stalled represented time returns no partial trajectory.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn simulate_viscous_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    areas_m2: &[f64],
    duration_s: f64,
    core: f64,
    law: PseDiffusion,
    options: ViscousAdaptiveOptions,
    budget: AdaptiveBudget,
) -> Result<AdaptiveRun, AdaptiveError> {
    let time = options.time;
    validate_options(duration_s, time)?;
    if !options.absolute_position_tolerance_m.is_finite()
        || options.absolute_position_tolerance_m <= 0.0
    {
        return Err(invalid_scalar(
            "adaptive absolute position tolerance",
            options.absolute_position_tolerance_m,
            None,
        )
        .into());
    }
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None).into());
    }
    let n = particles.len();
    diffusion::validate_inputs(n, areas_m2, duration_s)?;
    admit_cap("particles", n, budget.max_particles)?;
    let n64 = count(n, "adaptive moving particles")?;
    let active = duration_s > 0.0 && n >= 2;
    let diffuse = active && law.kinematic_viscosity_m2_s() > 0.0;
    let pairs_per_trial = if active {
        diffusion::pair_work(n, 1, if diffuse { 12 } else { 3 })?
    } else {
        0
    };
    let fractions_per_trial = if diffuse {
        diffusion::pair_work(n, 1, 7)?
    } else {
        0
    };
    let trial_work = n64
        .checked_mul(if diffuse { 4 } else { 2 })
        .and_then(|value| value.checked_add(pairs_per_trial))
        .and_then(|value| value.checked_add(if diffuse { 5 } else { 4 }))
        .ok_or(overflow("adaptive moving trial work"))?;
    admit_count("adaptive work", n64, budget.max_work_units)?;
    let peak_live_bytes = n
        .checked_mul(core::mem::size_of::<VortexParticle>())
        .and_then(|value| value.checked_mul(if active { 3 } else { 1 }))
        .ok_or(overflow("adaptive moving live bytes"))?;
    admit_cap("logical live bytes", peak_live_bytes, budget.max_live_bytes)?;
    let planned = if active { budget.max_work_units } else { n64 };
    let mut accountant = AdmittedBudget::admit_ambient(cx, planned)?;
    accountant.checkpoint("fs-vpm.viscous-adaptive-admission", cx)?;
    let mut meter = WorkMeter::default();
    let mut state = reserved_vec("viscous-adaptive-state", n)?;
    let mut coarse = reserved_vec("viscous-adaptive-coarse", if active { n } else { 0 })?;
    let mut fine = reserved_vec("viscous-adaptive-fine", if active { n } else { 0 })?;
    for (index, particle) in particles.iter().copied().enumerate() {
        meter.unit(&mut accountant, cx, "fs-vpm.viscous-adaptive-input", || {
            validate_particle(particle, index)?;
            diffusion::validate_area(areas_m2[index], index)?;
            state.push(particle);
            Ok(())
        })?;
    }
    meter.boundary(
        &mut accountant,
        cx,
        "fs-vpm.viscous-adaptive-input-complete",
    )?;

    let mut elapsed = if active { 0.0 } else { duration_s };
    let mut dt = time.initial_dt_s;
    let mut attempts = 0_usize;
    let mut accepted_steps = 0_usize;
    let mut rejected_steps = 0_usize;
    let mut pair_evaluations = 0_u64;
    let mut fraction_evaluations = 0_u64;
    let mut expected_work = n64;
    let mut max_accepted_error_ratio = 0.0_f64;
    while elapsed < duration_s {
        let next_attempt = attempts
            .checked_add(1)
            .ok_or(overflow("adaptive attempts"))?;
        admit_cap("adaptive attempts", next_attempt, budget.max_attempts)?;
        let next_pairs = pair_evaluations
            .checked_add(pairs_per_trial)
            .ok_or(overflow("adaptive moving pair work"))?;
        admit_count("pair evaluations", next_pairs, budget.max_pair_evaluations)?;
        let next_work = expected_work
            .checked_add(trial_work)
            .ok_or(overflow("adaptive moving accumulated work"))?;
        admit_count("adaptive work", next_work, budget.max_work_units)?;
        let (h, endpoint) = trial_time(elapsed, dt, duration_s)?;
        coarse.clear();
        fine.clear();
        for &particle in &state {
            meter.unit(&mut accountant, cx, "fs-vpm.viscous-adaptive-copy", || {
                coarse.push(particle);
                fine.push(particle);
                Ok(())
            })?;
        }
        // Keep h separate from its fractions so a small physical interval
        // survives multiplication by a large, representable original rate.
        super::viscous::weighted_step(
            &mut coarse,
            areas_m2,
            h,
            1.0,
            core,
            law,
            &mut pair_evaluations,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        super::viscous::weighted_step(
            &mut fine,
            areas_m2,
            h,
            0.5,
            core,
            law,
            &mut pair_evaluations,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        super::viscous::weighted_step(
            &mut fine,
            areas_m2,
            h,
            0.5,
            core,
            law,
            &mut pair_evaluations,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        let mut error_ratio = 0.0_f64;
        for i in 0..n {
            meter.unit(&mut accountant, cx, "fs-vpm.viscous-adaptive-error", || {
                let strength_error = scaled_error(
                    state[i].circulation,
                    coarse[i].circulation,
                    fine[i].circulation,
                    areas_m2[i],
                    time.absolute_vorticity_tolerance_s_inv,
                    time.relative_tolerance,
                );
                let position_error = position_error(
                    coarse[i].pos,
                    fine[i].pos,
                    options.absolute_position_tolerance_m,
                );
                error_ratio = error_ratio.max(strength_error).max(position_error);
                Ok(())
            })?;
        }
        if diffuse {
            // Step doubling alone can agree spuriously when both candidates
            // diffuse away the initial motion. Reverse the splitting order,
            // using the now-spare coarse buffer and the same physical owners.
            coarse.clear();
            for &particle in &state {
                meter.unit(&mut accountant, cx, "fs-vpm.viscous-witness-copy", || {
                    coarse.push(particle);
                    Ok(())
                })?;
            }
            super::viscous::reverse_step(
                &mut coarse,
                areas_m2,
                h,
                core,
                law,
                &mut pair_evaluations,
                &mut meter,
                &mut accountant,
                cx,
            )?;
            for i in 0..n {
                meter.unit(&mut accountant, cx, "fs-vpm.viscous-witness-error", || {
                    let strength_error = scaled_error(
                        state[i].circulation,
                        coarse[i].circulation,
                        fine[i].circulation,
                        areas_m2[i],
                        time.absolute_vorticity_tolerance_s_inv,
                        time.relative_tolerance,
                    );
                    let position_error = position_error(
                        coarse[i].pos,
                        fine[i].pos,
                        options.absolute_position_tolerance_m,
                    );
                    // This is an independent splitting witness, not a
                    // step-doubling estimate: undo the shared helpers' /3.
                    error_ratio = error_ratio
                        .max(3.0 * strength_error)
                        .max(3.0 * position_error);
                    Ok(())
                })?;
            }
        }
        meter.unit(&mut accountant, cx, "fs-vpm.viscous-adaptive-trial", || {
            Ok(())
        })?;
        meter.boundary(
            &mut accountant,
            cx,
            "fs-vpm.viscous-adaptive-trial-complete",
        )?;
        attempts = next_attempt;
        expected_work = next_work;
        fraction_evaluations = fraction_evaluations
            .checked_add(fractions_per_trial)
            .ok_or(overflow("adaptive moving kernel evaluations"))?;
        if pair_evaluations != next_pairs || meter.completed != expected_work {
            return Err(VpmError::AccountingMismatch {
                planned: expected_work,
                completed: meter.completed,
                charged: accountant.consumption().cost_charged,
            }
            .into());
        }
        if error_ratio <= 1.0 {
            core::mem::swap(&mut state, &mut fine);
            elapsed = endpoint;
            accepted_steps += 1;
            max_accepted_error_ratio = max_accepted_error_ratio.max(error_ratio);
        } else {
            rejected_steps += 1;
        }
        dt = next_step(time, h, error_ratio, elapsed)?;
    }
    meter.flush(&mut accountant, "fs-vpm.viscous-adaptive-finalize")?;
    if meter.completed != expected_work || accountant.consumption().cost_charged != expected_work {
        return Err(VpmError::AccountingMismatch {
            planned: expected_work,
            completed: meter.completed,
            charged: accountant.consumption().cost_charged,
        }
        .into());
    }
    accountant.checkpoint("fs-vpm.viscous-adaptive-publication", cx)?;
    Ok(AdaptiveRun {
        particles: state,
        elapsed_s: duration_s,
        accepted_steps,
        rejected_steps,
        attempts,
        pair_evaluations,
        fraction_evaluations,
        work_units: expected_work,
        peak_live_bytes,
        max_accepted_error_ratio,
        budget: accountant.consumption(),
    })
}

fn position_error(coarse: [f64; 2], fine: [f64; 2], tolerance: f64) -> f64 {
    let x = scaled_error(0.0, coarse[0], fine[0], 1.0, tolerance, 0.0);
    let y = scaled_error(0.0, coarse[1], fine[1], 1.0, tolerance, 0.0);
    let scale = x.max(y);
    if scale == 0.0 || !scale.is_finite() {
        return scale;
    }
    scale * det::sqrt((x / scale) * (x / scale) + (y / scale) * (y / scale))
}
