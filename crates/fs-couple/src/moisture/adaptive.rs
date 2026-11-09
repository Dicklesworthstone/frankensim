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
//! convergence/increment/stage-range failures and temporal error permit shorter
//! retries.

mod finite;
mod layered;
mod paired;
mod radial;
mod shared;
mod solution_accuracy;
mod solution_pair;

pub use solution_accuracy::{KineticPhaseAccuracy, SolutionMoistureTimeControls};
pub use solution_pair::{
    PairedSolutionMoistureAdvanceReport, PairedSolutionMoistureTimeControls,
    SolutionMoistureAdvanceAccount,
};

pub use layered::{
    LayeredMoistureAdvanceReport, MoistureSpatialPairedSolutionAdvanceReport,
    MoistureSpatialPairedSolutionTimeControls,
};
pub use paired::{AtmosphereMoistureAdvanceAccount, PairedAtmosphereMoistureAdvanceReport};

pub use shared::{
    SharedAtmosphereMoistureAdvanceCandidate, SharedCartesianAtmosphereMoistureAdvanceCandidate,
    SharedCartesianSolutionMoistureAdvanceCandidate, SharedMoistureAdvanceCandidate,
    SharedSolutionMoistureAdvanceCandidate,
};

use super::binding::MoistureField;
use super::cartesian::{MoistureCartesian, MoistureCartesianReport, MoistureCartesianSdirk2Report};
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
    /// failed/rejected probes. SDIRK2 charges two implicit stages per call,
    /// reserved before starting even if the call fails in its first stage.
    /// A complete comparison costs three ordinary calls or six SDIRK2 stages;
    /// each call retains its original scalar/nonlinear/linear work bounds.
    /// Must lie in 3..=100_000; SDIRK2 needs at least six to finish a comparison.
    pub maximum_evaluations: usize,
    pub moisture_absolute_tolerance: f64,
    /// Absolute error in dimensionless RH or solution water activity.
    pub activity_absolute_tolerance: f64,
    /// Resolves external water even when solution activity is nearly constant.
    /// Radial fields also use this for the L1 sum of shell water discrepancies.
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
    /// Accepted physical half intervals, excluding internal RK stages.
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    /// Charged endpoint/stage work, as in `maximum_evaluations`.
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

impl<Step, const N: usize> MoistureAdvanceReport<Step, N> {
    fn map_step<Next>(self, map: impl FnOnce(Step) -> Next) -> MoistureAdvanceReport<Next, N> {
        MoistureAdvanceReport {
            interval_s: self.interval_s,
            comparisons: self.comparisons,
            accepted_intervals: self.accepted_intervals,
            accepted_substeps: self.accepted_substeps,
            rejected_intervals: self.rejected_intervals,
            evaluations: self.evaluations,
            minimum_accepted_substep_s: self.minimum_accepted_substep_s,
            maximum_accepted_substep_s: self.maximum_accepted_substep_s,
            maximum_accepted_error_ratio: self.maximum_accepted_error_ratio,
            suggested_step_s: self.suggested_step_s,
            boundary_water_kg: self.boundary_water_kg,
            stored_water_kg: self.stored_water_kg,
            water_change_kg: self.water_change_kg,
            reservoir_water_kg: self.reservoir_water_kg,
            activity: self.activity,
            transport_balance_residual_kg: self.transport_balance_residual_kg,
            water_balance_residual_kg: self.water_balance_residual_kg,
            maximum_carrier_mass_residual_kg: self.maximum_carrier_mass_residual_kg,
            total_mass_residual_kg: self.total_mass_residual_kg,
            last_substep: map(self.last_substep),
        }
    }
}

/// A complete-window account for one original specimen, in its own face order.
#[derive(Clone, Debug)]
pub struct MoistureSpecimenAdvanceReport<const FACES: usize> {
    pub id: u64,
    pub boundary_water_kg: [f64; FACES],
    pub stored_water_kg: [f64; 2],
    pub water_change_kg: f64,
    pub transport_balance_residual_kg: f64,
    /// Maximum accepted componentwise moisture error, without other specimens
    /// or reservoir coordinates diluting this specimen's comparison.
    pub maximum_accepted_error_ratio: f64,
}

/// One shared-reservoir duration with independently accountable native fields.
#[derive(Clone, Debug)]
pub struct SharedMoistureAdvanceReport<Step, const FACES: usize> {
    /// Whole-system accounts; face totals sum the corresponding local face
    /// slots, not a reconstructed common surface or a fictitious geometry.
    pub summary: MoistureAdvanceReport<Step, FACES>,
    /// Stable-ID order. Each entry retains its own signed accepted transfers.
    pub specimens: Vec<MoistureSpecimenAdvanceReport<FACES>>,
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
    SpecimenWaterBalance {
        id: u64,
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
                | Self::StageRange { .. }
                | Self::ActivityStageRange { .. }
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
impl WaterStep<2> for MoistureSdirk2Report {
    fn boundary_water(&self) -> [f64; 2] {
        self.endpoint.boundary_water_kg
    }
}
impl WaterStep<6> for MoistureCartesianReport {
    fn boundary_water(&self) -> [f64; 6] {
        self.boundary_water_kg
    }
}
impl WaterStep<6> for MoistureCartesianSdirk2Report {
    fn boundary_water(&self) -> [f64; 6] {
        self.endpoint.boundary_water_kg
    }
}

#[derive(Clone, Copy)]
enum TimeMethod {
    FirstOrder,
    Sdirk2,
}
impl TimeMethod {
    fn evaluation_cost(self) -> usize {
        match self {
            Self::FirstOrder => 1,
            Self::Sdirk2 => 2,
        }
    }

    // For local order p+1, the fine error is the full/fine difference times
    // (u^(p+1)+v^(p+1))/(1-u^(p+1)-v^(p+1)), with u+v=1.
    fn fine_error_factor(self, u: f64, v: f64) -> f64 {
        match self {
            Self::FirstOrder => (u * u + v * v) / (2.0 * u * v),
            Self::Sdirk2 => (u * u * u + v * v * v) / (3.0 * u * v),
        }
    }

    fn doubling_threshold(self) -> f64 {
        // Doubling a span multiplies its leading local error by 2^(p+1).
        match self {
            Self::FirstOrder => 0.25,
            Self::Sdirk2 => 0.125,
        }
    }
}

// The single-field path keeps its view and transfer array inline. A shared
// collection holds one view per genuine field, never concatenated cell storage.
enum Specimens<T> {
    One(T),
    Many(Vec<T>),
}
impl<T> Specimens<T> {
    fn as_slice(&self) -> &[T] {
        match self {
            Self::One(value) => std::slice::from_ref(value),
            Self::Many(values) => values,
        }
    }
}
trait ReportTransfers<const N: usize> {
    fn transfers(&self) -> Specimens<(u64, [f64; N])>;
}
struct OneReport<Step>(Step);
impl<Step: WaterStep<N>, const N: usize> ReportTransfers<N> for OneReport<Step> {
    fn transfers(&self) -> Specimens<(u64, [f64; N])> {
        Specimens::One((0, self.0.boundary_water()))
    }
}

struct Snapshot<'a> {
    fields: Specimens<(u64, MoistureField<'a>)>,
    /// Optional material activity coordinates, in the same native field order.
    /// A shallow sorption curve must not hide a large chemical driving change
    /// behind a small dry-basis moisture increment.
    field_activities: Option<Specimens<(u64, &'a [f64])>>,
    reservoir_water_kg: Option<f64>,
    activity: Option<f64>,
}
fn field_snapshot<F: ReservoirField<N>, const N: usize>(
    field: &F,
) -> Result<Snapshot<'_>, MoistureError> {
    Ok(Snapshot {
        fields: Specimens::One((0, field.inventory())),
        field_activities: None,
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
fn snapshot_clock<E>(snapshot: &Snapshot<'_>, cx: &Cx<'_>) -> Result<f64, MoistureAdvanceError<E>> {
    let fields = snapshot.fields.as_slice();
    let Some((_, first)) = fields.first() else {
        return Err(MoistureAdvanceError::Invalid(
            "adaptive moisture needs at least one native field",
        ));
    };
    let time = first.time_s;
    if !time.is_finite() || fields.windows(2).any(|pair| pair[0].0 >= pair[1].0) {
        return Err(MoistureAdvanceError::Invalid(
            "adaptive specimen IDs must be unique and sorted, with a finite clock",
        ));
    }
    for (i, (_, field)) in fields.iter().enumerate() {
        if i % 128 == 0 {
            checkpoint(cx)?;
        }
        if field.time_s.to_bits() != time.to_bits() {
            return Err(MoistureAdvanceError::TimeResolution {
                time_s: field.time_s,
                target_s: time,
            });
        }
    }
    Ok(time)
}
fn check_clock<E>(
    snapshot: &Snapshot<'_>,
    time: f64,
    cx: &Cx<'_>,
) -> Result<(), MoistureAdvanceError<E>> {
    let actual = snapshot_clock(snapshot, cx)?;
    if actual.to_bits() != time.to_bits() {
        return Err(MoistureAdvanceError::TimeResolution {
            time_s: actual,
            target_s: time,
        });
    }
    Ok(())
}
fn same_field(a: &(u64, MoistureField<'_>), b: &(u64, MoistureField<'_>)) -> bool {
    a.0 == b.0
        && a.1.identity == b.1.identity
        && a.1.moisture.len() == b.1.moisture.len()
        && a.1.same_dry_masses(&b.1)
}
fn distance<E>(
    old: &Snapshot<'_>,
    a: &Snapshot<'_>,
    b: &Snapshot<'_>,
    c: MoistureTimeControls,
    field_errors: &mut [f64],
    cx: &Cx<'_>,
) -> Result<f64, MoistureAdvanceError<E>> {
    let activities = match (
        &old.field_activities,
        &a.field_activities,
        &b.field_activities,
    ) {
        (None, None, None) => None,
        (Some(old), Some(a), Some(b)) => Some((old.as_slice(), a.as_slice(), b.as_slice())),
        _ => {
            return Err(MoistureAdvanceError::Invalid(
                "material activity coordinates changed during one moisture duration",
            ));
        }
    };
    let (old, coarse, fine) = (
        old.fields.as_slice(),
        a.fields.as_slice(),
        b.fields.as_slice(),
    );
    if old.len() != coarse.len() || old.len() != fine.len() || old.len() != field_errors.len() {
        return Err(MoistureAdvanceError::Invalid(
            "adaptive comparison changed the specimen collection",
        ));
    }
    if activities.is_some_and(|(a, b, c)| {
        a.len() != old.len() || b.len() != old.len() || c.len() != old.len()
    }) {
        return Err(MoistureAdvanceError::Invalid(
            "material activity coordinates lost a native field",
        ));
    }
    let mut ratio: f64 = 0.0;
    for (slot, ((before, coarse), fine)) in old.iter().zip(coarse).zip(fine).enumerate() {
        checkpoint(cx)?;
        if !same_field(before, coarse) || !same_field(before, fine) {
            return Err(MoistureAdvanceError::Invalid(
                "adaptive comparison changed an original moisture geometry/source binding",
            ));
        }
        let mut local: f64 = 0.0;
        for (i, ((&before, &coarse), &fine)) in before
            .1
            .moisture
            .iter()
            .zip(coarse.1.moisture)
            .zip(fine.1.moisture)
            .enumerate()
        {
            if i % 128 == 0 {
                checkpoint(cx)?;
            }
            local = local.max(scaled_error(
                before,
                coarse,
                fine,
                c.moisture_absolute_tolerance,
                c.relative_tolerance,
            )?);
        }
        if let Some((old_activity, coarse_activity, fine_activity)) = activities {
            let coordinates = [
                &old_activity[slot],
                &coarse_activity[slot],
                &fine_activity[slot],
            ];
            if coordinates
                .iter()
                .any(|(id, values)| *id != before.0 || values.len() != before.1.moisture.len())
            {
                return Err(MoistureAdvanceError::Invalid(
                    "material activity coordinates changed their original cells",
                ));
            }
            for (i, ((&old, &a), &b)) in coordinates[0]
                .1
                .iter()
                .zip(coordinates[1].1)
                .zip(coordinates[2].1)
                .enumerate()
            {
                if i % 128 == 0 {
                    checkpoint(cx)?;
                }
                local = local.max(scaled_error(
                    old,
                    a,
                    b,
                    c.activity_absolute_tolerance,
                    c.relative_tolerance,
                )?);
            }
        }
        field_errors[slot] = local;
        ratio = ratio.max(local);
    }
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
    cost: usize,
    propose: &mut impl FnMut(&S, f64) -> Result<(S, R), E>,
    cx: &Cx<'_>,
) -> Result<(S, R), MoistureAdvanceError<E>> {
    checkpoint(cx)?;
    if cost > maximum - *calls {
        return Err(MoistureAdvanceError::WorkLimit {
            time_s: time,
            accepted_intervals: accepted,
            evaluations: *calls,
        });
    }
    *calls += cost;
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
    let (state, report) = evolve_fields(
        initial,
        target_s,
        controls,
        |state, dt| propose(state, dt).map(|(state, step)| (state, OneReport(step))),
        snapshot,
        cx,
    )?;
    Ok((state, report.summary.map_step(|step| step.0)))
}

fn evolve_with_error<S: Clone, R: WaterStep<N>, E: Retryable, const N: usize>(
    initial: &S,
    target_s: f64,
    controls: MoistureTimeControls,
    mut propose: impl FnMut(&S, f64) -> Result<(S, R), E>,
    snapshot: impl for<'a> Fn(&'a S) -> Result<Snapshot<'a>, E>,
    additional_error: impl Fn(&S, &S, &S) -> Result<f64, MoistureAdvanceError<E>>,
    cx: &Cx<'_>,
) -> Result<(S, MoistureAdvanceReport<R, N>), MoistureAdvanceError<E>> {
    let (state, report) = evolve_fields_with_error(
        initial,
        target_s,
        controls,
        |state, dt| propose(state, dt).map(|(state, step)| (state, OneReport(step))),
        snapshot,
        additional_error,
        cx,
    )?;
    Ok((state, report.summary.map_step(|step| step.0)))
}

fn evolve_fields<S: Clone, R: ReportTransfers<N>, E: Retryable, const N: usize>(
    initial: &S,
    target_s: f64,
    controls: MoistureTimeControls,
    propose: impl FnMut(&S, f64) -> Result<(S, R), E>,
    snapshot: impl for<'a> Fn(&'a S) -> Result<Snapshot<'a>, E>,
    cx: &Cx<'_>,
) -> Result<(S, SharedMoistureAdvanceReport<R, N>), MoistureAdvanceError<E>> {
    evolve_fields_with_error(
        initial,
        target_s,
        controls,
        propose,
        snapshot,
        |_, _, _| Ok(0.0),
        cx,
    )
}

// Paired chambers need separate RH and water coordinates: their opposite
// exchanges can cancel in the total reservoir inventory and specimen storage.
// They still use this one native full-versus-two-half controller and accounts.
#[allow(clippy::too_many_arguments)]
fn evolve_fields_with_error<S: Clone, R: ReportTransfers<N>, E: Retryable, const N: usize>(
    initial: &S,
    target_s: f64,
    controls: MoistureTimeControls,
    propose: impl FnMut(&S, f64) -> Result<(S, R), E>,
    snapshot: impl for<'a> Fn(&'a S) -> Result<Snapshot<'a>, E>,
    additional_error: impl Fn(&S, &S, &S) -> Result<f64, MoistureAdvanceError<E>>,
    cx: &Cx<'_>,
) -> Result<(S, SharedMoistureAdvanceReport<R, N>), MoistureAdvanceError<E>> {
    evolve_fields_with_method(
        initial,
        target_s,
        controls,
        TimeMethod::FirstOrder,
        propose,
        snapshot,
        additional_error,
        cx,
    )
}

#[allow(clippy::too_many_arguments)]
fn evolve_fields_with_method<S: Clone, R: ReportTransfers<N>, E: Retryable, const N: usize>(
    initial: &S,
    target_s: f64,
    controls: MoistureTimeControls,
    method: TimeMethod,
    mut propose: impl FnMut(&S, f64) -> Result<(S, R), E>,
    snapshot: impl for<'a> Fn(&'a S) -> Result<Snapshot<'a>, E>,
    additional_error: impl Fn(&S, &S, &S) -> Result<f64, MoistureAdvanceError<E>>,
    cx: &Cx<'_>,
) -> Result<(S, SharedMoistureAdvanceReport<R, N>), MoistureAdvanceError<E>> {
    checkpoint(cx)?;
    let original = snapshot(initial).map_err(MoistureAdvanceError::FinalInventory)?;
    let start = snapshot_clock(&original, cx)?;
    validate(start, target_s, controls)?;
    let mut state = initial.clone();
    let mut time = start;
    let mut next_dt = controls.initial_step_s;
    let (mut comparisons, mut accepted, mut rejected, mut evaluations) = (0, 0, 0, 0);
    let (mut minimum, mut maximum, mut maximum_error) = (f64::INFINITY, 0.0_f64, 0.0_f64);
    let mut water = vec![[WaterSum::default(); N]; original.fields.as_slice().len()];
    let mut field_errors = vec![0.0; water.len()];
    let mut field_maximum_error = vec![0.0_f64; water.len()];
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
                method.evaluation_cost(),
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
                method.evaluation_cost(),
                &mut propose,
                cx,
            )?;
            check_clock(
                &snapshot(&half).map_err(MoistureAdvanceError::FinalInventory)?,
                middle,
                cx,
            )?;
            let first_water = first.transfers();
            drop(first);
            let (fine, report) = probe(
                &half,
                middle,
                h2,
                accepted,
                &mut evaluations,
                controls.maximum_evaluations,
                method.evaluation_cost(),
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
        check_clock(&a, end, cx)?;
        check_clock(&b, end, cx)?;
        // Method-order local error for possibly unequal representable halves.
        let (u, v) = (h1 / span, h2 / span);
        let factor = method.fine_error_factor(u, v);
        let field_error = distance(&before, &a, &b, controls, &mut field_errors, cx)?;
        let activity_error = optional_error(
            before.activity,
            a.activity,
            b.activity,
            controls.activity_absolute_tolerance,
            controls.relative_tolerance,
        )?;
        let reservoir_error = optional_error(
            before.reservoir_water_kg,
            a.reservoir_water_kg,
            b.reservoir_water_kg,
            controls.reservoir_water_absolute_tolerance_kg,
            controls.relative_tolerance,
        )?;
        let extra = additional_error(&state, &coarse, &fine)?;
        if !extra.is_finite() || extra < 0.0 {
            return Err(MoistureAdvanceError::Invalid(
                "unrepresentable additional moisture temporal comparison",
            ));
        }
        let error = factor
            * field_error
                .max(activity_error)
                .max(reservoir_error)
                .max(extra);
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
        let second_water = report.transfers();
        let (first_water, second_water) = (first_water.as_slice(), second_water.as_slice());
        if first_water.len() != water.len() || second_water.len() != water.len() {
            return Err(MoistureAdvanceError::Invalid(
                "native moisture receipt lost a specimen's face transfers",
            ));
        }
        for (i, sums) in water.iter_mut().enumerate() {
            checkpoint(cx)?;
            let id = original.fields.as_slice()[i].0;
            if first_water[i].0 != id || second_water[i].0 != id {
                return Err(MoistureAdvanceError::Invalid(
                    "native moisture receipt changed stable specimen IDs",
                ));
            }
            for (side, sum) in sums.iter_mut().enumerate() {
                sum.add(first_water[i].1[side]);
                sum.add(second_water[i].1[side]);
            }
            field_maximum_error[i] = field_maximum_error[i].max(factor * field_errors[i]);
        }
        accepted += 1;
        minimum = minimum.min(h1.min(h2));
        maximum = maximum.max(h1.max(h2));
        maximum_error = maximum_error.max(error);
        last = Some(report);
        state = fine;
        time = end;
        next_dt = if error <= method.doubling_threshold() {
            span + span.min(controls.maximum_step_s - span)
        } else {
            span
        };
        next_dt = next_dt.clamp(controls.minimum_step_s, controls.maximum_step_s);
    }
    let final_state = snapshot(&state).map_err(MoistureAdvanceError::FinalInventory)?;
    check_clock(&final_state, target_s, cx)?;
    let (original_fields, final_fields) =
        (original.fields.as_slice(), final_state.fields.as_slice());
    if original_fields.len() != final_fields.len() {
        return Err(MoistureAdvanceError::Invalid(
            "final moisture endpoint lost a specimen",
        ));
    }
    let mut specimens = Vec::with_capacity(water.len());
    let mut changed_water = WaterSum::default();
    let mut stored = [WaterSum::default(); 2];
    let mut all_faces = [WaterSum::default(); N];
    for (slot, (old, new)) in original_fields.iter().zip(final_fields).enumerate() {
        checkpoint(cx)?;
        if !same_field(old, new) {
            return Err(MoistureAdvanceError::Invalid(
                "final moisture endpoint changed an original specimen binding",
            ));
        }
        let mut changed = WaterSum::default();
        for (i, (&after, &before)) in new.1.moisture.iter().zip(old.1.moisture).enumerate() {
            if i % 128 == 0 {
                checkpoint(cx)?;
            }
            changed.add((after - before) * old.1.dry_mass_at(i));
        }
        let change = changed.total();
        let transfers = water[slot].map(WaterSum::total);
        let mut transferred = WaterSum::default();
        for (sum, value) in all_faces.iter_mut().zip(transfers) {
            sum.add(value);
            transferred.add(value);
        }
        let residual = change - transferred.total();
        let stored_water_kg = [old.1.water_mass_kg(), new.1.water_mass_kg()];
        if stored_water_kg.iter().any(|m| !m.is_finite()) {
            return Err(MoistureAdvanceError::Invalid(
                "unrepresentable original/final specimen water inventory",
            ));
        }
        if original_fields.len() > 1
            && (!residual.is_finite() || residual.abs() > controls.water_balance_tolerance_kg)
        {
            return Err(MoistureAdvanceError::SpecimenWaterBalance {
                id: old.0,
                residual_kg: residual,
                tolerance_kg: controls.water_balance_tolerance_kg,
            });
        }
        changed_water.add(change);
        for (sum, value) in stored.iter_mut().zip(stored_water_kg) {
            sum.add(value);
        }
        specimens.push(MoistureSpecimenAdvanceReport {
            id: old.0,
            boundary_water_kg: transfers,
            stored_water_kg,
            water_change_kg: change,
            transport_balance_residual_kg: residual,
            maximum_accepted_error_ratio: field_maximum_error[slot],
        });
    }
    let change = changed_water.total();
    let transfers = all_faces.map(WaterSum::total);
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
    let stored_water = stored.map(WaterSum::total);
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
    Ok((
        state,
        SharedMoistureAdvanceReport {
            summary: report,
            specimens,
        },
    ))
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

    /// Evolve a complete duration with the original second-order SDIRK2 stages.
    /// Faces and source temperature stay fixed. Only the two-half endpoint and
    /// its RK-weighted water transfers are accepted; failed prefixes stay private.
    ///
    /// Equal halves estimate the fine local error as one third of the full/fine
    /// difference. The dyadic controller uses cubic local-error scaling, with
    /// the corresponding factor for unequal representable halves. This is an
    /// asymptotic local estimate, not a global error bound. Physical `StageRange`
    /// refusals shorten the span without clipping or falling back to first order.
    ///
    /// Every attempted SDIRK2 call charges two stages against the work limit,
    /// including early failures; a complete comparison costs six. Reports count
    /// accepted physical half intervals separately from those implicit stages.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, MoistureAdvanceReport<MoistureSdirk2Report, 2>),
        MoistureAdvanceError<MoistureError>,
    > {
        let (state, report) = evolve_fields_with_method(
            self,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                state
                    .advance_sdirk2(dt, faces, step, nonlinear, cx)
                    .map(|(state, report)| (state, OneReport(report)))
            },
            field_snapshot::<Self, 2>,
            |_, _, _| Ok(0.0),
            cx,
        )?;
        Ok((state, report.summary.map_step(|step| step.0)))
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

    /// Evolve a complete 3D duration with second-order implicit source stages.
    ///
    /// The original full-versus-two-half controller uses cubic local-error
    /// scaling and the correct factor for unequal representable halves.
    /// Each attempted endpoint charges two implicit stages, including failed
    /// attempts. Only accepted physical half intervals and their RK-weighted
    /// six-face transfers enter the whole-duration water account.
    ///
    /// Physical stage-range failures shorten the proposal without clipping or
    /// a first-order fallback. A later failure publishes no accepted prefix.
    /// Temperature and external faces remain fixed over the complete call;
    /// this is an asymptotic local error estimate, not a global error bound
    /// or an adaptive finite-reservoir RK coupling.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureFace; 6],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureAdvanceReport<MoistureCartesianSdirk2Report, 6>,
        ),
        MoistureAdvanceError<MoistureError>,
    > {
        let (state, report) = evolve_fields_with_method(
            self,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                state
                    .advance_sdirk2(dt, faces, step, nonlinear, cx)
                    .map(|(state, report)| (state, OneReport(report)))
            },
            field_snapshot::<Self, 6>,
            |_, _, _| Ok(0.0),
            cx,
        )?;
        Ok((state, report.summary.map_step(|step| step.0)))
    }
}
