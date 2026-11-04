//! Adaptive positive-time windows of the original reversible vortex maps.
//!
//! Compare one step with two half steps and retain the fine state without
//! extrapolation. Every accepted pair subflow therefore keeps its original
//! circulation and linear/angular impulse invariants. A held affine source
//! retains its own actual moment account; rejected candidates add no physical
//! source moment. The local estimate controls temporal splitting only.

use super::{AdaptiveBudget, AdaptiveError, AdaptiveRun};
use crate::conservative::ambient::{AffineFlow2, AmbientMomentChange, HalfFlow};
use crate::conservative::{Composition, CompositionOrder, weighted_symmetric_sweep};
use crate::diffusion::adaptive::{admit_count, count, next_step_with_growth, overflow};
use crate::{
    AdmittedBudget, Cx, VortexParticle, VpmError, WorkMeter, admit_cap, invalid_scalar,
    reserved_vec, validate_particle,
};
use fs_math::{dd::Dd, det};

/// Physical step limits and a local position-error target.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveOptions {
    /// First positive trial, between the declared minimum and maximum [s].
    pub initial_dt_s: f64,
    /// Smallest retry [s]. The final endpoint remainder may be shorter.
    pub min_dt_s: f64,
    /// Largest physical trial [s]. Fourth-order internal stages may be signed.
    pub max_dt_s: f64,
    /// Positive absolute local position tolerance [m].
    pub absolute_position_tolerance_m: f64,
    /// Nonnegative relative tolerance times the larger initial/fine cloud
    /// radius about its unweighted centroid. No world-origin magnitude or
    /// division by net circulation enters the scale.
    pub relative_tolerance: f64,
}

/// One whole-window execution budget plus the explicit affine-source cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AmbientAdaptiveBudget {
    /// Original adaptive particle/attempt/pair/work/live-payload envelope.
    pub run: AdaptiveBudget,
    /// All source visits, including rejected and coarse candidates.
    pub max_source_particle_updates: u64,
}

/// Accepted adaptive window and the physical source increments on its path.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientAdaptiveRun {
    pub run: AdaptiveRun,
    /// Only accepted fine trajectories contribute to this original DD account.
    pub ambient_moment_change: AmbientMomentChange,
    /// Actual computational source visits, including rejected/coarse trials.
    pub source_particle_updates: u64,
}

/// Complete an inviscid window with controlled local position error.
///
/// For order p=2 or 4, acceptance requires
/// `max_i |x_fine-x_coarse| / ((2^p-1)*(atol+rtol*max(R_old,R_fine))) <= 1`.
/// R is the maximum Euclidean radius about the unweighted cloud centroid.
/// Both the error norm and its length scale are invariant under rigid changes
/// of world frame up to floating-point representation. Zero net circulation
/// is admitted. The original positive core/noncolliding physical assumptions
/// remain required; close interactions must still be spatially resolved.
///
/// Rejected trials halve dt. Accepted ratios <= 1/8 (Second) or 1/32 (Fourth)
/// permit doubling. Physical endpoint time uses compensated accumulation;
/// fractional half steps remain separate from dt, retaining the original
/// scaled maps even when dt is subnormal. No Richardson state correction,
/// moment projection, viscosity, remapping or source adaptation is introduced.
/// This is a local splitting-error estimate, not a certified/global trajectory
/// bound. Variable accepted timesteps do not establish symplecticity or time
/// reversibility, and conservation alone does not establish phase accuracy.
///
/// Three particle banks are reserved once and reused. Initial validation and
/// centered geometry cost 2 N units; each active trial costs 3 q N(N-1) pair
/// visits plus 3 N+4 units, with q=1 or 3 for Second/Fourth. Empty, zero-time
/// and unforced singleton windows are identities with one bank and N units.
/// Rejected trials count against the original whole-window work/attempt caps.
///
/// # Errors
/// Invalid physical/options/resource data, an original numerical refusal,
/// cancellation, exhausted aggregate caps, a tolerance floor, unrepresentable
/// centered geometry or stalled time. No partial trajectory is returned.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    duration_s: f64,
    core: f64,
    order: CompositionOrder,
    options: AdaptiveOptions,
    budget: AdaptiveBudget,
) -> Result<AdaptiveRun, AdaptiveError> {
    run(
        cx,
        particles,
        duration_s,
        core,
        order,
        None,
        options,
        AmbientAdaptiveBudget {
            run: budget,
            max_source_particle_updates: 0,
        },
    )
    .map(|result| result.run)
}

/// The same adaptive window with one held incompressible affine source.
///
/// All source maps are the original exact `HalfFlow` owner. Adjacent maps at
/// the fine trajectory's midpoint coalesce: a trial uses 5 N source visits at
/// Second order and 11 N at Fourth, rather than 6 N and 12 N. The full and
/// half weighted maps are prepared once per trial on the stack. No trial
/// allocates. Work includes all actual source visits; physical moment changes
/// include only accepted fine motion. Uniform flow does not inflate the
/// relative tolerance through a large world position or common translation.
///
/// The finite-volume/body momentum and continuum-error limitations of the
/// original prescribed-source owner remain. A zero field takes the exact
/// unforced path; no diffusion law can enter the signed internal composition.
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_affine_flow_cx(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    duration_s: f64,
    core: f64,
    order: CompositionOrder,
    flow: AffineFlow2,
    options: AdaptiveOptions,
    budget: AmbientAdaptiveBudget,
) -> Result<AmbientAdaptiveRun, AdaptiveError> {
    run(
        cx,
        particles,
        duration_s,
        core,
        order,
        (!flow.is_zero()).then_some(flow),
        options,
        budget,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run(
    cx: &Cx<'_>,
    particles: &[VortexParticle],
    duration_s: f64,
    core: f64,
    order: CompositionOrder,
    flow: Option<AffineFlow2>,
    options: AdaptiveOptions,
    budget: AmbientAdaptiveBudget,
) -> Result<AmbientAdaptiveRun, AdaptiveError> {
    validate_options(duration_s, core, options)?;
    let n = particles.len();
    admit_cap("particles", n, budget.run.max_particles)?;
    let active = duration_s > 0.0 && n > 0 && (n > 1 || flow.is_some());
    let n64 = count(n, "adaptive particles")?;
    let initial_work = n64
        .checked_mul(if active { 2 } else { 1 })
        .ok_or(overflow("adaptive input work"))?;
    let pairs_per_trial = if active {
        crate::diffusion::pair_work(n, 3, order.sweeps() as u64)?
    } else {
        0
    };
    let sources_per_trial = if active && flow.is_some() {
        n64.checked_mul(3 * order.source_passes() - 1)
            .ok_or(overflow("adaptive source visits"))?
    } else {
        0
    };
    let trial_work = n64
        .checked_mul(3)
        .and_then(|v| v.checked_add(4))
        .and_then(|v| v.checked_add(pairs_per_trial))
        .and_then(|v| v.checked_add(sources_per_trial))
        .ok_or(overflow("adaptive trial work"))?;
    admit_count("adaptive work", initial_work, budget.run.max_work_units)?;
    let peak_live_bytes = n
        .checked_mul(core::mem::size_of::<VortexParticle>())
        .and_then(|v| v.checked_mul(if active { 3 } else { 1 }))
        .ok_or(overflow("adaptive live bytes"))?;
    admit_cap(
        "logical live bytes",
        peak_live_bytes,
        budget.run.max_live_bytes,
    )?;
    let mut accountant = AdmittedBudget::admit_ambient(
        cx,
        if active {
            budget.run.max_work_units
        } else {
            initial_work
        },
    )?;
    accountant.checkpoint("fs-vpm.adaptive-admission", cx)?;
    let mut state = reserved_vec("adaptive-state", n)?;
    let mut coarse = reserved_vec("adaptive-coarse", if active { n } else { 0 })?;
    let mut fine = reserved_vec("adaptive-fine", if active { n } else { 0 })?;
    let mut meter = WorkMeter::default();
    let mut center = [Dd::ZERO; 2];
    let divisor = dd(n.max(1) as f64);
    for (i, particle) in particles.iter().copied().enumerate() {
        meter.unit(&mut accountant, cx, "fs-vpm.adaptive-input", || {
            validate_particle(particle, i)?;
            if active {
                add_center(&mut center, particle, divisor)?;
            }
            state.push(particle);
            Ok(())
        })?;
    }
    let mut radius = 0.0_f64;
    if active {
        for (i, particle) in state.iter().enumerate() {
            meter.unit(&mut accountant, cx, "fs-vpm.adaptive-geometry", || {
                radius = radius.max(centered_radius(*particle, center, i)?);
                Ok(())
            })?;
        }
    }
    meter.boundary(&mut accountant, cx, "fs-vpm.adaptive-input-complete")?;
    let composition = Composition::new(order);
    let (error_divisor, growth_ratio) = match order {
        CompositionOrder::Second => (3.0, 0.125),
        CompositionOrder::Fourth => (15.0, 0.03125),
    };
    let end = dd(duration_s);
    let mut elapsed = if active { Dd::ZERO } else { end };
    let mut dt = options.initial_dt_s;
    let (mut attempts, mut accepted, mut rejected) = (0_usize, 0_usize, 0_usize);
    let (mut pairs, mut source_updates) = (0_u64, 0_u64);
    let mut expected_work = initial_work;
    let mut moments = [Dd::ZERO; 3];
    let mut max_accepted_error_ratio = 0.0_f64;
    while positive(end - elapsed) {
        let next_attempt = attempts
            .checked_add(1)
            .ok_or(overflow("adaptive attempts"))?;
        admit_cap("adaptive attempts", next_attempt, budget.run.max_attempts)?;
        let next_pairs = pairs
            .checked_add(pairs_per_trial)
            .ok_or(overflow("adaptive accumulated pairs"))?;
        admit_count(
            "pair evaluations",
            next_pairs,
            budget.run.max_pair_evaluations,
        )?;
        let next_sources = source_updates
            .checked_add(sources_per_trial)
            .ok_or(overflow("adaptive accumulated sources"))?;
        admit_count(
            "source particle updates",
            next_sources,
            budget.max_source_particle_updates,
        )?;
        let next_work = expected_work
            .checked_add(trial_work)
            .ok_or(overflow("adaptive accumulated work"))?;
        admit_count("adaptive work", next_work, budget.run.max_work_units)?;
        let (h, endpoint) = compensated_trial_time(elapsed, dt, end)?;
        // Fallible map preparation precedes all candidate work. The previous
        // input/trial boundary has flushed the meter, so an unrepresentable
        // source map cannot leave completed copy visits uncharged.
        let full_source = flow
            .map(|f| SourceMaps::new(f, h, &composition, 1.0))
            .transpose()?;
        let half_source = flow
            .map(|f| SourceMaps::new(f, h, &composition, 0.5))
            .transpose()?;
        coarse.clear();
        fine.clear();
        for &particle in &state {
            meter.unit(&mut accountant, cx, "fs-vpm.adaptive-copy", || {
                coarse.push(particle);
                fine.push(particle);
                Ok(())
            })?;
        }
        let mut coarse_moments = [Dd::ZERO; 3];
        let mut fine_moments = [Dd::ZERO; 3];
        evolve(
            &mut coarse,
            h,
            core,
            &composition,
            1,
            full_source.as_ref(),
            None,
            &mut pairs,
            &mut source_updates,
            &mut coarse_moments,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        evolve(
            &mut fine,
            h,
            core,
            &composition,
            2,
            half_source.as_ref(),
            full_source.as_ref().map(|s| &s.outer),
            &mut pairs,
            &mut source_updates,
            &mut fine_moments,
            &mut meter,
            &mut accountant,
            cx,
        )?;
        let mut fine_center = [Dd::ZERO; 2];
        for &particle in &fine {
            meter.unit(&mut accountant, cx, "fs-vpm.adaptive-center", || {
                add_center(&mut fine_center, particle, divisor)
            })?;
        }
        let (mut fine_radius, mut max_error) = (0.0_f64, 0.0_f64);
        for i in 0..n {
            meter.unit(&mut accountant, cx, "fs-vpm.adaptive-error", || {
                fine_radius = fine_radius.max(centered_radius(fine[i], fine_center, i)?);
                let delta = core::array::from_fn::<_, 2, _>(|d| {
                    (dd(fine[i].pos[d]) - dd(coarse[i].pos[d])).to_f64()
                });
                let error = det::hypot(delta[0], delta[1]);
                if !error.is_finite() {
                    return Err(range_error("fs-vpm.adaptive-error", i));
                }
                max_error = max_error.max(error);
                Ok(())
            })?;
        }
        let error_ratio =
            normalized_error(max_error, radius.max(fine_radius), options, error_divisor);
        meter.unit(&mut accountant, cx, "fs-vpm.adaptive-trial", || Ok(()))?;
        meter.boundary(&mut accountant, cx, "fs-vpm.adaptive-trial-complete")?;
        if error_ratio.is_nan() || pairs != next_pairs || source_updates != next_sources {
            return Err(range_error("fs-vpm.adaptive-trial", 0).into());
        }
        attempts = next_attempt;
        expected_work = next_work;
        if error_ratio <= 1.0 {
            for d in 0..3 {
                moments[d] = moments[d] + fine_moments[d];
                if !finite(moments[d]) {
                    return Err(range_error("fs-vpm.adaptive-source-sum", 0).into());
                }
            }
            core::mem::swap(&mut state, &mut fine);
            radius = fine_radius;
            elapsed = endpoint;
            accepted += 1;
            max_accepted_error_ratio = max_accepted_error_ratio.max(error_ratio);
        } else {
            rejected += 1;
        }
        dt = next_step_with_growth(
            options.min_dt_s,
            options.max_dt_s,
            growth_ratio,
            h,
            error_ratio,
            elapsed.to_f64(),
        )?;
    }
    meter.flush(&mut accountant, "fs-vpm.adaptive-finalize")?;
    let consumption = accountant.consumption();
    if meter.completed != expected_work || consumption.cost_charged != expected_work {
        return Err(VpmError::AccountingMismatch {
            planned: expected_work,
            completed: meter.completed,
            charged: consumption.cost_charged,
        }
        .into());
    }
    accountant.checkpoint("fs-vpm.adaptive-publication", cx)?;
    Ok(AmbientAdaptiveRun {
        run: AdaptiveRun {
            particles: state,
            elapsed_s: duration_s,
            accepted_steps: accepted,
            rejected_steps: rejected,
            attempts,
            pair_evaluations: pairs,
            fraction_evaluations: 0,
            work_units: expected_work,
            peak_live_bytes,
            max_accepted_error_ratio,
            budget: consumption,
        },
        ambient_moment_change: AmbientMomentChange {
            first_moment_m3_s: [moments[0].to_f64(), moments[1].to_f64()],
            quadratic_moment_m4_s: moments[2].to_f64(),
        },
        source_particle_updates: source_updates,
    })
}

struct SourceMaps {
    outer: HalfFlow,
    bridge: Option<HalfFlow>,
}
impl SourceMaps {
    fn new(
        flow: AffineFlow2,
        h: f64,
        composition: &Composition,
        fraction: f64,
    ) -> Result<Self, VpmError> {
        Ok(Self {
            outer: HalfFlow::weighted(flow, h, composition.weights()[0] * fraction)?,
            bridge: composition
                .bridge_weight()
                .map(|w| HalfFlow::weighted(flow, h, w * fraction))
                .transpose()?,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn evolve(
    state: &mut [VortexParticle],
    h: f64,
    core: f64,
    composition: &Composition,
    substeps: usize,
    source: Option<&SourceMaps>,
    join: Option<&HalfFlow>,
    pairs: &mut u64,
    source_updates: &mut u64,
    moments: &mut [Dd; 3],
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    if let Some(source) = source {
        apply(
            &source.outer,
            state,
            source_updates,
            moments,
            meter,
            accountant,
            cx,
        )?;
    }
    let fraction = if substeps == 2 { 0.5 } else { 1.0 };
    for step in 0..substeps {
        for (i, &weight) in composition.weights().iter().enumerate() {
            if i > 0
                && let Some(bridge) = source.and_then(|s| s.bridge.as_ref())
            {
                apply(
                    bridge,
                    state,
                    source_updates,
                    moments,
                    meter,
                    accountant,
                    cx,
                )?;
            }
            weighted_symmetric_sweep(
                state,
                h,
                core,
                weight * fraction,
                pairs,
                meter,
                accountant,
                cx,
            )?;
        }
        meter.unit(accountant, cx, "fs-vpm.conservative-step", || Ok(()))?;
        if step + 1 < substeps
            && let Some(join) = join
        {
            apply(join, state, source_updates, moments, meter, accountant, cx)?;
        }
    }
    if let Some(source) = source {
        apply(
            &source.outer,
            state,
            source_updates,
            moments,
            meter,
            accountant,
            cx,
        )?;
    }
    Ok(())
}

fn apply(
    source: &HalfFlow,
    state: &mut [VortexParticle],
    updates: &mut u64,
    moments: &mut [Dd; 3],
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
) -> Result<(), VpmError> {
    source.apply(state, moments, meter, accountant, cx)?;
    *updates += state.len() as u64; // The whole trial's checked plan admits this sum.
    Ok(())
}

fn add_center(center: &mut [Dd; 2], particle: VortexParticle, divisor: Dd) -> Result<(), VpmError> {
    for d in 0..2 {
        center[d] = center[d] + dd(particle.pos[d]) / divisor;
        if !finite(center[d]) {
            return Err(range_error("fs-vpm.adaptive-center", 0));
        }
    }
    Ok(())
}

fn centered_radius(p: VortexParticle, center: [Dd; 2], i: usize) -> Result<f64, VpmError> {
    let d = core::array::from_fn::<_, 2, _>(|a| (dd(p.pos[a]) - center[a]).to_f64());
    let r = det::hypot(d[0], d[1]);
    if !r.is_finite() {
        return Err(range_error("fs-vpm.adaptive-radius", i));
    }
    Ok(r)
}

fn normalized_error(error: f64, radius: f64, options: AdaptiveOptions, divisor: f64) -> f64 {
    if error == 0.0 {
        return 0.0;
    }
    let absolute = options.absolute_position_tolerance_m;
    let scale = error.max(radius).max(absolute);
    let numerator = error / scale;
    let denominator = absolute / scale + options.relative_tolerance * (radius / scale);
    let ratio = (numerator / denominator) / divisor;
    if numerator > 0.0 && denominator.is_normal() && ratio.is_finite() {
        return ratio;
    }
    // Preserve tiny absolute tolerances even beside a very large common cloud
    // scale; an underflowed normalized 0/0 must not admit an inaccurate trial.
    let a = det::ln(absolute);
    let b = if options.relative_tolerance > 0.0 && radius > 0.0 {
        det::ln(options.relative_tolerance) + det::ln(radius)
    } else {
        f64::NEG_INFINITY
    };
    let largest = a.max(b);
    let log_denominator = largest + det::ln(det::exp(a - largest) + det::exp(b - largest));
    det::exp(det::ln(error) - det::ln(divisor) - log_denominator)
}

fn validate_options(duration: f64, core: f64, o: AdaptiveOptions) -> Result<(), VpmError> {
    if !duration.is_finite() || duration < 0.0 {
        return Err(invalid_scalar("adaptive duration", duration, None));
    }
    if !core.is_finite() || core < 0.0 {
        return Err(invalid_scalar("core", core, None));
    }
    for (field, value) in [
        ("adaptive initial dt", o.initial_dt_s),
        ("adaptive minimum dt", o.min_dt_s),
        ("adaptive maximum dt", o.max_dt_s),
        (
            "adaptive absolute position tolerance",
            o.absolute_position_tolerance_m,
        ),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid_scalar(field, value, None));
        }
    }
    if o.min_dt_s > o.initial_dt_s || o.initial_dt_s > o.max_dt_s {
        return Err(invalid_scalar("adaptive dt ordering", o.initial_dt_s, None));
    }
    if !o.relative_tolerance.is_finite() || o.relative_tolerance < 0.0 {
        return Err(invalid_scalar(
            "adaptive relative tolerance",
            o.relative_tolerance,
            None,
        ));
    }
    Ok(())
}

// Shared by source-coupled positive viscous windows. Retain every low time
// part instead of inventing a tiny final source interval from repeated f64
// addition. The original inviscid arithmetic and refusal boundary are intact.
pub(super) fn compensated_trial_time(
    elapsed: Dd,
    dt: f64,
    end: Dd,
) -> Result<(f64, Dd), AdaptiveError> {
    let remaining = end - elapsed;
    let mut h = dt.min(remaining.to_f64());
    // A rounded-up endpoint remainder must not advance beyond the window.
    // Its remaining low part is retained as a final physical interval.
    if positive(dd(h) - remaining) && h > 0.0 {
        h = f64::from_bits(h.to_bits() - 1);
    }
    let endpoint = elapsed + dd(h);
    if h <= 0.0 || !finite(endpoint) || !positive(endpoint - elapsed) || positive(endpoint - end) {
        return Err(AdaptiveError::TimeStagnation {
            elapsed_bits: elapsed.to_f64().to_bits(),
            dt_bits: h.to_bits(),
        });
    }
    Ok((h, endpoint))
}

fn dd(v: f64) -> Dd {
    Dd::from_f64(v)
}
fn finite(v: Dd) -> bool {
    v.hi.is_finite() && v.lo.is_finite() && v.to_f64().is_finite()
}
pub(super) fn positive(v: Dd) -> bool {
    v.hi > 0.0 || (v.hi == 0.0 && v.lo > 0.0)
}
fn range_error(phase: &'static str, particle: usize) -> VpmError {
    VpmError::NonFiniteResult { phase, particle }
}
