//! Whole-duration accuracy around the original conservative moisture endpoints.
//!
//! One full interval estimates the error of two half intervals. Only those
//! original half steps contribute water or accepted state; no extrapolated
//! moisture, reservoir or source state is constructed. Componentwise content,
//! external water and RH/activity discrepancies use explicit tolerances on
//! changes. The estimate is local and asymptotic, not a global error bound.
//!
//! Boundaries and source temperature stay fixed over one call; callers cut
//! intervals at environmental changes. Every failed call discards the entire
//! proposed duration, including successful internal prefixes. Only typed
//! convergence/increment failures and temporal error permit shorter retries.

mod finite;

use super::binding::MoistureField;
use super::cartesian::{MoistureCartesian, MoistureCartesianReport};
use super::reservoir::ReservoirField;
use super::*;

/// Temporal accuracy and whole-call work/inventory limits, independent of each
/// original endpoint's linear, nonlinear and physical conservation controls.
#[derive(Clone, Copy, Debug)]
pub struct MoistureTimeControls {
    pub initial_step_s: f64,
    /// Floor for comparison spans, not their accepted physical half steps.
    /// A shorter final remainder may be attempted once, without bypassing error.
    pub minimum_step_s: f64,
    /// Largest proposed duration; actual intervals use representable endpoints.
    pub maximum_step_s: f64,
    /// Counts every native field or coupled-reservoir endpoint call, including
    /// failed/rejected probes. Each comparison uses at most three; each native
    /// call retains its original scalar/nonlinear/linear work bounds.
    /// Must lie in 3..=100_000.
    pub maximum_evaluations: usize,
    pub moisture_absolute_tolerance: f64,
    /// Absolute error in dimensionless RH or solution water activity.
    pub activity_absolute_tolerance: f64,
    /// Resolves external water even when solution activity is nearly constant.
    pub reservoir_water_absolute_tolerance_kg: f64,
    /// Multiplies each coordinate's change from the last accepted state; [0,1).
    pub relative_tolerance: f64,
    /// One fixed original-to-final water budget [kg], never multiplied by steps.
    /// Finite vapor also checks whole-window total and per-carrier mass here.
    pub water_balance_tolerance_kg: f64,
}

/// Only accepted fine transfers contribute to the whole-duration accounts.
/// Source history is bounded: the last native half-step receipt is retained.
#[derive(Clone, Debug)]
pub struct MoistureAdvanceReport<Step, const FACES: usize> {
    pub interval_s: [f64; 2],
    pub comparisons: usize,
    pub accepted_intervals: usize,
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    pub maximum_accepted_error_ratio: f64,
    pub suggested_step_s: f64,
    pub boundary_water_kg: [f64; FACES],
    pub stored_water_kg: [f64; 2],
    /// Compensated sum of original-to-final local dry-mass water changes.
    pub water_change_kg: f64,
    pub reservoir_water_kg: Option<[f64; 2]>,
    pub activity: Option<[f64; 2]>,
    /// Original-to-final bound water minus integrated accepted boundary water.
    pub transport_balance_residual_kg: f64,
    /// Includes external reservoir water when present; otherwise equals above.
    pub water_balance_residual_kg: f64,
    /// Present for finite vapor; independently rechecked over the whole duration.
    pub maximum_carrier_mass_residual_kg: Option<f64>,
    pub total_mass_residual_kg: Option<f64>,
    pub last_substep: Step,
}

/// No variant exposes a partially advanced physical field or reservoir.
#[derive(Debug)]
pub enum MoistureAdvanceError<E> {
    Invalid(&'static str),
    Step {
        time_s: f64,
        step_s: f64,
        evaluations: usize,
        error: E,
    },
    FinalInventory(E),
    Accuracy {
        time_s: f64,
        step_s: f64,
        error_ratio: f64,
    },
    WorkLimit {
        time_s: f64,
        accepted_intervals: usize,
        evaluations: usize,
    },
    TimeResolution {
        time_s: f64,
        target_s: f64,
    },
    WaterBalance {
        residual_kg: f64,
        tolerance_kg: f64,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for MoistureAdvanceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive moisture: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MoistureAdvanceError<E> {}

trait Retryable {
    fn retryable(&self) -> bool;
}
impl Retryable for MoistureError {
    fn retryable(&self) -> bool {
        matches!(
            self,
            Self::NotConverged(_)
                | Self::NonlinearNotConverged(_)
                | Self::SpatialNotConverged { .. }
                | Self::ContentIncrement { .. }
        )
    }
}
trait WaterStep<const N: usize> {
    fn boundary_water(&self) -> [f64; N];
}
impl WaterStep<2> for MoistureStepReport {
    fn boundary_water(&self) -> [f64; 2] {
        self.boundary_water_kg
    }
}
impl WaterStep<6> for MoistureCartesianReport {
    fn boundary_water(&self) -> [f64; 6] {
        self.boundary_water_kg
    }
}

struct Snapshot<'a> {
    field: MoistureField<'a>,
    reservoir_water_kg: Option<f64>,
    activity: Option<f64>,
}
fn field_snapshot<F: ReservoirField<N>, const N: usize>(
    field: &F,
) -> Result<Snapshot<'_>, MoistureError> {
    Ok(Snapshot {
        field: field.inventory(),
        reservoir_water_kg: None,
        activity: None,
    })
}

// Compensated accepted-face totals; rejected probes never enter this account.
#[derive(Clone, Copy, Default)]
struct WaterSum {
    sum: f64,
    correction: f64,
}
impl WaterSum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn total(self) -> f64 {
        self.sum + self.correction
    }
}

fn checkpoint<E>(cx: &Cx<'_>) -> Result<(), MoistureAdvanceError<E>> {
    cx.checkpoint().map_err(|_| MoistureAdvanceError::Cancelled)
}
fn validate<E>(
    start: f64,
    target: f64,
    c: MoistureTimeControls,
) -> Result<(), MoistureAdvanceError<E>> {
    if !target.is_finite()
        || target <= start
        || !positive(target - start)
        || !positive(c.minimum_step_s)
        || !positive(c.maximum_step_s)
        || c.minimum_step_s > c.maximum_step_s
        || !positive(c.initial_step_s)
        || c.initial_step_s < c.minimum_step_s
        || c.initial_step_s > c.maximum_step_s
        || !(3..=100_000).contains(&c.maximum_evaluations)
        || !positive(c.moisture_absolute_tolerance)
        || !positive(c.activity_absolute_tolerance)
        || !positive(c.reservoir_water_absolute_tolerance_kg)
        || !c.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.relative_tolerance)
        || !positive(c.water_balance_tolerance_kg)
    {
        return Err(MoistureAdvanceError::Invalid(
            "forward target, positive temporal/inventory scales and bounded timestep/work controls required",
        ));
    }
    Ok(())
}
fn scaled_error<E>(
    old: f64,
    a: f64,
    b: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, MoistureAdvanceError<E>> {
    let allowance = absolute + relative * (a - old).abs().max((b - old).abs());
    if !positive(allowance) || !a.is_finite() || !b.is_finite() || !old.is_finite() {
        return Err(MoistureAdvanceError::Invalid(
            "unrepresentable moisture temporal comparison",
        ));
    }
    Ok((a - b).abs() / allowance)
}
fn optional_error<E>(
    old: Option<f64>,
    a: Option<f64>,
    b: Option<f64>,
    absolute: f64,
    relative: f64,
) -> Result<f64, MoistureAdvanceError<E>> {
    match (old, a, b) {
        (None, None, None) => Ok(0.0),
        (Some(old), Some(a), Some(b)) => scaled_error(old, a, b, absolute, relative),
        _ => Err(MoistureAdvanceError::Invalid(
            "reservoir observable changed during one moisture duration",
        )),
    }
}
fn distance<E>(
    old: &Snapshot<'_>,
    a: &Snapshot<'_>,
    b: &Snapshot<'_>,
    c: MoistureTimeControls,
    cx: &Cx<'_>,
) -> Result<f64, MoistureAdvanceError<E>> {
    if old.field.identity != a.field.identity
        || old.field.identity != b.field.identity
        || old.field.moisture.len() != a.field.moisture.len()
        || old.field.moisture.len() != b.field.moisture.len()
    {
        return Err(MoistureAdvanceError::Invalid(
            "adaptive comparison changed the original moisture geometry/source binding",
        ));
    }
    let mut ratio: f64 = 0.0;
    for (i, ((&before, &coarse), &fine)) in old
        .field
        .moisture
        .iter()
        .zip(a.field.moisture)
        .zip(b.field.moisture)
        .enumerate()
    {
        if i % 128 == 0 {
            checkpoint(cx)?;
        }
        ratio = ratio.max(scaled_error(
            before,
            coarse,
            fine,
            c.moisture_absolute_tolerance,
            c.relative_tolerance,
        )?);
    }
    ratio = ratio.max(optional_error(
        old.activity,
        a.activity,
        b.activity,
        c.activity_absolute_tolerance,
        c.relative_tolerance,
    )?);
    ratio = ratio.max(optional_error(
        old.reservoir_water_kg,
        a.reservoir_water_kg,
        b.reservoir_water_kg,
        c.reservoir_water_absolute_tolerance_kg,
        c.relative_tolerance,
    )?);
    Ok(ratio)
}

#[allow(clippy::too_many_arguments)]
fn probe<S, R, E>(
    state: &S,
    time: f64,
    dt: f64,
    accepted: usize,
    calls: &mut usize,
    maximum: usize,
    propose: &mut impl FnMut(&S, f64) -> Result<(S, R), E>,
    cx: &Cx<'_>,
) -> Result<(S, R), MoistureAdvanceError<E>> {
    checkpoint(cx)?;
    if *calls == maximum {
        return Err(MoistureAdvanceError::WorkLimit {
            time_s: time,
            accepted_intervals: accepted,
            evaluations: *calls,
        });
    }
    *calls += 1;
    let result = propose(state, dt);
    checkpoint(cx)?;
    result.map_err(|error| MoistureAdvanceError::Step {
        time_s: time,
        step_s: dt,
        evaluations: *calls,
        error,
    })
}

fn evolve<S: Clone, R: WaterStep<N>, E: Retryable, const N: usize>(
    initial: &S,
    target_s: f64,
    controls: MoistureTimeControls,
    mut propose: impl FnMut(&S, f64) -> Result<(S, R), E>,
    snapshot: impl for<'a> Fn(&'a S) -> Result<Snapshot<'a>, E>,
    cx: &Cx<'_>,
) -> Result<(S, MoistureAdvanceReport<R, N>), MoistureAdvanceError<E>> {
    checkpoint(cx)?;
    let original = snapshot(initial).map_err(MoistureAdvanceError::FinalInventory)?;
    let start = original.field.time_s;
    validate(start, target_s, controls)?;
    let mut state = initial.clone();
    let mut time = start;
    let mut next_dt = controls.initial_step_s;
    let (mut comparisons, mut accepted, mut rejected, mut evaluations) = (0, 0, 0, 0);
    let (mut minimum, mut maximum, mut maximum_error) = (f64::INFINITY, 0.0_f64, 0.0_f64);
    let mut water = [WaterSum::default(); N];
    let mut last = None;
    while time < target_s {
        checkpoint(cx)?;
        let end = if next_dt >= target_s - time {
            target_s
        } else {
            time + next_dt
        };
        let span = end - time;
        let middle = time + 0.5 * span;
        if !(time < middle && middle < end) {
            return Err(MoistureAdvanceError::TimeResolution {
                time_s: time,
                target_s: end,
            });
        }
        let (h1, h2) = (middle - time, end - middle);
        comparisons += 1;
        let pair = (|| {
            let (coarse, _) = probe(
                &state,
                time,
                span,
                accepted,
                &mut evaluations,
                controls.maximum_evaluations,
                &mut propose,
                cx,
            )?;
            let (half, first) = probe(
                &state,
                time,
                h1,
                accepted,
                &mut evaluations,
                controls.maximum_evaluations,
                &mut propose,
                cx,
            )?;
            if snapshot(&half)
                .map_err(MoistureAdvanceError::FinalInventory)?
                .field
                .time_s
                .to_bits()
                != middle.to_bits()
            {
                return Err(MoistureAdvanceError::TimeResolution {
                    time_s: time,
                    target_s: middle,
                });
            }
            let first_water = first.boundary_water();
            drop(first);
            let (fine, report) = probe(
                &half,
                middle,
                h2,
                accepted,
                &mut evaluations,
                controls.maximum_evaluations,
                &mut propose,
                cx,
            )?;
            Ok((coarse, fine, report, first_water))
        })();
        let (coarse, fine, report, first_water) = match pair {
            Ok(result) => result,
            Err(error) => {
                let retry =
                    matches!(&error, MoistureAdvanceError::Step { error, .. } if error.retryable());
                let smaller = (span * 0.5).max(controls.minimum_step_s);
                if retry && smaller < span {
                    next_dt = smaller;
                    rejected += 1;
                    continue;
                }
                return Err(error);
            }
        };
        let before = snapshot(&state).map_err(MoistureAdvanceError::FinalInventory)?;
        let a = snapshot(&coarse).map_err(MoistureAdvanceError::FinalInventory)?;
        let b = snapshot(&fine).map_err(MoistureAdvanceError::FinalInventory)?;
        if a.field.time_s.to_bits() != end.to_bits() || b.field.time_s.to_bits() != end.to_bits() {
            return Err(MoistureAdvanceError::TimeResolution {
                time_s: time,
                target_s: end,
            });
        }
        // First-order local error for possibly unequal representable halves.
        let (u, v) = (h1 / span, h2 / span);
        let factor = (u * u + v * v) / (2.0 * u * v);
        let error = factor * distance(&before, &a, &b, controls, cx)?;
        if error > 1.0 {
            let smaller = (span * 0.5).max(controls.minimum_step_s);
            if smaller < span {
                next_dt = smaller;
                rejected += 1;
                continue;
            }
            return Err(MoistureAdvanceError::Accuracy {
                time_s: time,
                step_s: span,
                error_ratio: error,
            });
        }
        checkpoint(cx)?;
        for (sum, (a, b)) in water
            .iter_mut()
            .zip(first_water.into_iter().zip(report.boundary_water()))
        {
            sum.add(a);
            sum.add(b);
        }
        accepted += 1;
        minimum = minimum.min(h1.min(h2));
        maximum = maximum.max(h1.max(h2));
        maximum_error = maximum_error.max(error);
        last = Some(report);
        state = fine;
        time = end;
        next_dt = if error <= 0.25 {
            span + span.min(controls.maximum_step_s - span)
        } else {
            span
        };
        next_dt = next_dt.clamp(controls.minimum_step_s, controls.maximum_step_s);
    }
    let final_state = snapshot(&state).map_err(MoistureAdvanceError::FinalInventory)?;
    let mut changed_water = WaterSum::default();
    for (i, (new, old)) in final_state
        .field
        .moisture
        .iter()
        .zip(original.field.moisture)
        .enumerate()
    {
        if i % 128 == 0 {
            checkpoint(cx)?;
        }
        changed_water.add((new - old) * original.field.dry_mass_per_cell_kg);
    }
    let change = changed_water.total();
    let transfers = water.map(WaterSum::total);
    let mut transferred_water = WaterSum::default();
    for value in transfers {
        transferred_water.add(value);
    }
    let transport_balance = change - transferred_water.total();
    let reservoir_water = original
        .reservoir_water_kg
        .zip(final_state.reservoir_water_kg)
        .map(|(a, b)| [a, b]);
    let balance = reservoir_water.map_or(transport_balance, |[a, b]| b - a + change);
    let stored_water = [
        original.field.water_mass_kg(),
        final_state.field.water_mass_kg(),
    ];
    if stored_water
        .iter()
        .chain(reservoir_water.as_ref().into_iter().flatten())
        .any(|m| !m.is_finite())
    {
        return Err(MoistureAdvanceError::Invalid(
            "unrepresentable original/final moisture water inventory",
        ));
    }
    for residual in [transport_balance, balance] {
        if !residual.is_finite() || residual.abs() > controls.water_balance_tolerance_kg {
            return Err(MoistureAdvanceError::WaterBalance {
                residual_kg: residual,
                tolerance_kg: controls.water_balance_tolerance_kg,
            });
        }
    }
    let report = MoistureAdvanceReport {
        interval_s: [start, time],
        comparisons,
        accepted_intervals: accepted,
        accepted_substeps: 2 * accepted,
        rejected_intervals: rejected,
        evaluations,
        minimum_accepted_substep_s: minimum,
        maximum_accepted_substep_s: maximum,
        maximum_accepted_error_ratio: maximum_error,
        suggested_step_s: next_dt,
        boundary_water_kg: transfers,
        stored_water_kg: stored_water,
        water_change_kg: change,
        reservoir_water_kg: reservoir_water,
        activity: original
            .activity
            .zip(final_state.activity)
            .map(|(a, b)| [a, b]),
        transport_balance_residual_kg: transport_balance,
        water_balance_residual_kg: balance,
        maximum_carrier_mass_residual_kg: None,
        total_mass_residual_kg: None,
        last_substep: last.expect("positive completed duration contains accepted half steps"),
    };
    checkpoint(cx)?;
    Ok((state, report))
}

impl MoistureSlab {
    /// Evolve the complete positive duration with constant faces and the slab's
    /// configured lagged/implicit source coefficients. Only the final field is
    /// returned; accepted internal prefixes remain private on any later failure.
    pub fn advance_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureFace; 2],
        step: MoistureStepConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, MoistureAdvanceReport<MoistureStepReport, 2>),
        MoistureAdvanceError<MoistureError>,
    > {
        evolve(
            self,
            target_s,
            time,
            |state, dt| state.advance(dt, faces, step, cx),
            field_snapshot::<Self, 2>,
            cx,
        )
    }
}
impl MoistureCartesian {
    /// Complete 3D diffusion with componentwise temporal accuracy and a fixed
    /// whole-duration water budget. Source temperature and six faces stay fixed.
    pub fn advance_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureFace; 6],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, MoistureAdvanceReport<MoistureCartesianReport, 6>),
        MoistureAdvanceError<MoistureError>,
    > {
        evolve(
            self,
            target_s,
            time,
            |state, dt| state.advance(dt, faces, step, nonlinear, cx),
            field_snapshot::<Self, 6>,
            cx,
        )
    }
}
