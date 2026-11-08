//! Held electrical drive histories on the original implicit thermal owner.
//!
//! `TimeSignal` owns the SI schedule and right-continuous hold convention.
//! Each incoming interval retains its held drive through its endpoint; a jump
//! changes the next interval's drive without resetting the body's enthalpy.
//! The explicit maximum step and optional native step-doubling estimate control
//! temporal resolution independently of source-feedback admission and energy
//! residual tolerances. Accepted states always come from the original solver.

use super::{IncandescentError, IncandescentStep, PreparedIncandescentSource};
use fs_conduction::ConductionError;
use fs_conduction::lumped::MAX_LUMPED_ENTHALPY_STEPS;
use fs_conduction::lumped::endpoint::LumpedEndpointEnvironment;
use fs_material::conductor::{OhmicDrive, graph::MaterialOhmicNode};
use fs_material::phase::UniformEnthalpyStepInput;
use fs_qty::Dims;
use fs_scenario::{Interp, ScenarioError, TimeSignal};

/// Bind one existing SI signal to a held ideal current or voltage drive.
/// Only constants and tables with `Interp::Hold` are admitted. Linear ramps,
/// carrier waveforms and cycle-averaged powers are not inferred from samples.
#[derive(Clone, Copy, Debug)]
pub enum IncandescentDriveHistory<'a> {
    Current(&'a TimeSignal),
    Voltage(&'a TimeSignal),
}

impl<'a> IncandescentDriveHistory<'a> {
    fn signal(self) -> &'a TimeSignal {
        match self {
            Self::Current(signal) | Self::Voltage(signal) => signal,
        }
    }

    fn drive(self, value: f64) -> OhmicDrive {
        match self {
            Self::Current(_) => OhmicDrive::Current(value),
            Self::Voltage(_) => OhmicDrive::Voltage(value),
        }
    }

    fn dimensions(self) -> Dims {
        match self {
            Self::Current(_) => Dims([0, 0, 0, 0, 1, 0]),
            Self::Voltage(_) => Dims([2, 1, -3, 0, -1, 0]),
        }
    }
}

/// Physical clock and bounded thermal work, shared by both history modes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncandescentHistoryControls {
    /// Absolute start in the original signal's seconds-based clock.
    pub start_s: f64,
    /// Positive maximum physical substep [s], also clipped at every drive knot.
    pub maximum_step_s: f64,
    /// Bounds ALL native solves, including discarded coarse or rejected trials
    /// in adaptive mode, and retained rows. A held table may contain at most
    /// this many entries plus one, including entries outside the horizon.
    pub maximum_steps: usize,
}

/// One accepted interval with the original complete electrical/thermal report.
/// The electrical state uses the interval's incoming held drive, including
/// when `end_s` is exactly a drive discontinuity.
#[derive(Clone, Debug, PartialEq)]
pub struct TimedIncandescentStep {
    pub start_s: f64,
    pub end_s: f64,
    pub step: IncandescentStep,
}

/// Local time-error target in the original primary enthalpy coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncandescentHistoryAccuracy {
    /// Absolute specific-enthalpy difference between a full native step and
    /// two half steps [J/kg]. This detects phase evolution at fixed temperature
    /// and is independent of the chosen enthalpy datum. It must exceed eight
    /// times the combined native enthalpy tolerance and state roundoff scale.
    pub maximum_local_enthalpy_error_j_kg: f64,
}

/// Accepted native half steps and the work used to resolve their time history.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveIncandescentHistory {
    pub steps: Vec<TimedIncandescentStep>,
    /// Includes discarded coarse steps and all three solves of rejected trials.
    pub native_solves: usize,
    pub rejected_intervals: usize,
    /// Largest accepted full-versus-two-half-step enthalpy estimate [J/kg].
    /// This is a local estimate, not a bound on accumulated trajectory error.
    pub maximum_accepted_error_j_kg: f64,
}

#[derive(Debug)]
pub enum IncandescentHistoryError {
    Native(IncandescentError),
    Signal(ScenarioError),
    Invalid(&'static str),
}

impl core::fmt::Display for IncandescentHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => write!(f, "incandescent history: {error}"),
            Self::Signal(error) => write!(f, "incandescent drive history: {error}"),
            Self::Invalid(what) => write!(f, "incandescent history: {what}"),
        }
    }
}
impl std::error::Error for IncandescentHistoryError {}
impl From<IncandescentError> for IncandescentHistoryError {
    fn from(error: IncandescentError) -> Self {
        Self::Native(error)
    }
}
impl From<ConductionError> for IncandescentHistoryError {
    fn from(error: ConductionError) -> Self {
        Self::Native(error.into())
    }
}
impl From<ScenarioError> for IncandescentHistoryError {
    fn from(error: ScenarioError) -> Self {
        Self::Signal(error)
    }
}

/// Advance a complete held-drive history without changing the accepted input.
///
/// `input.duration_s` is the requested horizon. The represented clock horizon
/// must agree within `64 * f64::EPSILON` relative to that duration; a large
/// absolute origin cannot silently stretch or shrink it. `input.internal_heat_j`
/// supplies additional constant prescribed power, `heat / represented_duration`,
/// independently of Joule heating. Each represented subinterval receives that
/// power times its actual duration. The original source, geometry, caloric chart
/// and reservoir remain fixed. Phase-dependent electrical sources use the
/// original chart-aware preparation; no phase law is reconstructed here.
///
/// Signal validation occurs once through its native cancellable evaluator.
/// Subsequent traversal reads its immutable held entries directly, with no
/// interpolation, schedule clone or repeated whole-table validation. Native
/// feedback bounds are prepared once per held span, then every substep calls
/// `PreparedIncandescentSource::advance`. Time increments that cannot be
/// represented at the supplied clock origin refuse explicitly.
///
/// All rows publish together after the horizon completes. Cancellation, source
/// refusal, allocation failure or exhausted work returns no partial vector.
/// Energy closure remains native and per step; callers must refine the temporal
/// maximum step for their own temperature/light-history accuracy requirements.
pub fn advance_history(
    source: &MaterialOhmicNode,
    environment: &impl LumpedEndpointEnvironment,
    input: UniformEnthalpyStepInput<'_>,
    drive_history: IncandescentDriveHistory<'_>,
    controls: IncandescentHistoryControls,
    checkpoint: impl FnMut() -> Result<(), ConductionError>,
) -> Result<Vec<TimedIncandescentStep>, IncandescentHistoryError> {
    Ok(advance_history_impl(
        source,
        environment,
        input,
        drive_history,
        controls,
        None,
        checkpoint,
    )?
    .steps)
}

/// Resolve held-drive transients using the original implicit step's time error.
///
/// Each trial compares one full backward-Euler step with two half steps from
/// the same accepted state. For this first-order method the local fine-state
/// error estimate is `abs(h_fine - h_coarse) / (2^1 - 1)`. If clock rounding
/// makes the two durations `a, b` unequal, the leading-order fine-state estimate
/// is multiplied by `(a*a + b*b) / (2*a*b)`, evaluated without squaring times.
/// A midpoint must be strictly inside the represented interval. Accepted trials
/// retain BOTH native half-step reports; no extrapolated state or synthetic
/// source/boundary energy is published. Rejected trials halve the interval;
/// an estimate below one quarter of the target permits doubling the next one.
/// The physical maximum, native feedback bound and exact drive knots still
/// clip every trial. Native source refusals remain errors, not retry signals.
///
/// The absolute J/kg target remains meaningful through melting and freezing.
/// It must exceed the algebraic/roundoff scale; it is never inflated to make a
/// trial pass. Step doubling supplies a LOCAL estimate, not a rigorous global
/// temperature, phase, or emitted-light error bound, especially at caloric
/// corners. The clock, source and atomic-publication contracts of
/// [`advance_history`] also apply. `maximum_steps` counts every native solve.
pub fn advance_history_adaptive(
    source: &MaterialOhmicNode,
    environment: &impl LumpedEndpointEnvironment,
    input: UniformEnthalpyStepInput<'_>,
    drive_history: IncandescentDriveHistory<'_>,
    controls: IncandescentHistoryControls,
    accuracy: IncandescentHistoryAccuracy,
    checkpoint: impl FnMut() -> Result<(), ConductionError>,
) -> Result<AdaptiveIncandescentHistory, IncandescentHistoryError> {
    advance_history_impl(
        source,
        environment,
        input,
        drive_history,
        controls,
        Some(accuracy),
        checkpoint,
    )
}

fn advance_history_impl(
    source: &MaterialOhmicNode,
    environment: &impl LumpedEndpointEnvironment,
    input: UniformEnthalpyStepInput<'_>,
    drive_history: IncandescentDriveHistory<'_>,
    controls: IncandescentHistoryControls,
    accuracy: Option<IncandescentHistoryAccuracy>,
    mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
) -> Result<AdaptiveIncandescentHistory, IncandescentHistoryError> {
    use IncandescentHistoryError::Invalid;
    checkpoint()?;
    if !controls.start_s.is_finite()
        || !controls.maximum_step_s.is_finite()
        || controls.maximum_step_s <= 0.0
        || controls.maximum_steps == 0
        || controls.maximum_steps > MAX_LUMPED_ENTHALPY_STEPS
        || !input.duration_s.is_finite()
        || input.duration_s <= 0.0
        || !input.internal_heat_j.is_finite()
    {
        return Err(Invalid(
            "finite clock, positive horizon/step and bounded nonzero step budget required",
        ));
    }
    if let Some(accuracy) = accuracy
        && (!accuracy.maximum_local_enthalpy_error_j_kg.is_finite()
            || accuracy.maximum_local_enthalpy_error_j_kg <= 0.0
            || accuracy.maximum_local_enthalpy_error_j_kg
                <= 8.0 * environment.enthalpy_tolerance_j_kg())
    {
        return Err(Invalid(
            "local time-error target must be finite, positive and above eight native enthalpy tolerances",
        ));
    }
    let end = controls.start_s + input.duration_s;
    let represented_duration = end - controls.start_s;
    let prescribed_power = input.internal_heat_j / represented_duration;
    if !end.is_finite()
        || !represented_duration.is_finite()
        || represented_duration <= 0.0
        || (represented_duration / input.duration_s - 1.0).abs() > 64.0 * f64::EPSILON
        || !prescribed_power.is_finite()
        || (input.internal_heat_j != 0.0 && prescribed_power == 0.0)
    {
        return Err(Invalid(
            "physical horizon or prescribed source power is not representable at the requested precision",
        ));
    }
    let signal = drive_history.signal();
    let (times, values): (&[f64], &[f64]) = match signal {
        TimeSignal::Constant(_) => (&[], &[]),
        TimeSignal::Table {
            times,
            values,
            interp: Interp::Hold,
            ..
        } => (times, values),
        _ => {
            return Err(Invalid(
                "electrical history requires a constant or previous-sample hold table",
            ));
        }
    };
    if signal.dims() != drive_history.dimensions() {
        return Err(Invalid(
            "drive history must use amperes for current or volts for voltage",
        ));
    }
    if times.len() > controls.maximum_steps + 1 || values.len() > controls.maximum_steps + 1 {
        return Err(Invalid("held drive table exceeds the declared work budget"));
    }
    let mut held = signal
        .eval_with_checkpoint(controls.start_s, &mut |_| checkpoint())??
        .value;
    let mut next_event = times.partition_point(|&time| time <= controls.start_s);
    let mut time = controls.start_s;
    let mut current = input;
    let mut report = AdaptiveIncandescentHistory {
        steps: Vec::new(),
        native_solves: 0,
        rejected_intervals: 0,
        maximum_accepted_error_j_kg: 0.0,
    };
    while time < end {
        checkpoint()?;
        if report.native_solves == controls.maximum_steps {
            return Err(Invalid(
                "incandescent history exhausted its native step budget",
            ));
        }
        let span_end = times.get(next_event).copied().unwrap_or(end).min(end);
        let prepared = PreparedIncandescentSource::try_new_with_caloric_curve(
            source,
            drive_history.drive(held),
            input.curve,
        )?;
        let mut proposal = current;
        proposal.duration_s = controls.maximum_step_s.min(span_end - time);
        proposal.internal_heat_j = prescribed_power * proposal.duration_s;
        let limit = prepared
            .maximum_step_s(environment, proposal, &mut checkpoint)?
            .min(controls.maximum_step_s);
        let mut next_duration = limit;
        while time < span_end {
            checkpoint()?;
            let required_solves = if accuracy.is_some() { 3 } else { 1 };
            if controls.maximum_steps - report.native_solves < required_solves {
                return Err(Invalid(
                    "incandescent history exhausted its native step budget",
                ));
            }
            let remaining = span_end - time;
            let maximum_duration = next_duration.min(remaining);
            let mut next = if maximum_duration == remaining {
                span_end
            } else {
                time + maximum_duration
            };
            // Adding a small step to an absolute clock can round upward.
            // The duration passed to the native owner must still obey its bound.
            if next - time > maximum_duration {
                next = next.next_down();
            }
            let duration = next - time;
            if !next.is_finite() || !duration.is_finite() || duration <= 0.0 {
                return Err(Invalid(
                    "native substep cannot be represented on this clock",
                ));
            }
            report
                .steps
                .try_reserve(if accuracy.is_some() { 2 } else { 1 })
                .map_err(|_| Invalid("history output allocation failed"))?;
            if let Some(accuracy) = accuracy {
                let middle = time + 0.5 * duration;
                let first_duration = middle - time;
                let second_duration = next - middle;
                if middle <= time || middle >= next {
                    return Err(Invalid(
                        "adaptive midpoint cannot be represented inside this clock interval",
                    ));
                }
                let split_correction = 0.5 * (first_duration / second_duration)
                    + 0.5 * (second_duration / first_duration);
                if !split_correction.is_finite() || split_correction <= 0.0 {
                    return Err(Invalid("adaptive half-interval ratio is not representable"));
                }
                let split_correction = split_correction.max(1.0);
                let coarse = native_interval(
                    &prepared,
                    environment,
                    current,
                    duration,
                    prescribed_power,
                    &mut checkpoint,
                )?;
                let first = native_interval(
                    &prepared,
                    environment,
                    current,
                    first_duration,
                    prescribed_power,
                    &mut checkpoint,
                )?;
                let mut half_input = current;
                half_input.initial = first.thermal.state;
                let second = native_interval(
                    &prepared,
                    environment,
                    half_input,
                    second_duration,
                    prescribed_power,
                    &mut checkpoint,
                )?;
                report.native_solves += 3;
                let coarse_h = coarse.thermal.state.specific_enthalpy_j_kg();
                let fine_h = second.thermal.state.specific_enthalpy_j_kg();
                let state_scale = current
                    .initial
                    .specific_enthalpy_j_kg()
                    .abs()
                    .max(coarse_h.abs())
                    .max(fine_h.abs());
                let noise =
                    8.0 * (environment.enthalpy_tolerance_j_kg() + f64::EPSILON * state_scale);
                let tolerance = accuracy.maximum_local_enthalpy_error_j_kg;
                if tolerance <= noise {
                    return Err(Invalid(
                        "local time-error target is below combined native enthalpy and roundoff scales",
                    ));
                }
                let error = (fine_h - coarse_h).abs() * split_correction;
                if !error.is_finite() {
                    return Err(Invalid("local enthalpy time-error estimate is not finite"));
                }
                if error > tolerance {
                    report.rejected_intervals += 1;
                    next_duration = 0.5 * duration;
                    continue;
                }
                report.maximum_accepted_error_j_kg = report.maximum_accepted_error_j_kg.max(error);
                current.initial = second.thermal.state;
                report.steps.push(TimedIncandescentStep {
                    start_s: time,
                    end_s: middle,
                    step: first,
                });
                report.steps.push(TimedIncandescentStep {
                    start_s: middle,
                    end_s: next,
                    step: second,
                });
                next_duration = if error <= 0.25 * tolerance {
                    if duration >= 0.5 * limit {
                        limit
                    } else {
                        2.0 * duration
                    }
                } else {
                    duration
                };
            } else {
                let step = native_interval(
                    &prepared,
                    environment,
                    current,
                    duration,
                    prescribed_power,
                    &mut checkpoint,
                )?;
                report.native_solves += 1;
                current.initial = step.thermal.state;
                report.steps.push(TimedIncandescentStep {
                    start_s: time,
                    end_s: next,
                    step,
                });
            }
            time = next;
        }
        if time < end {
            held = values[next_event];
            next_event += 1;
        }
    }
    checkpoint()?;
    Ok(report)
}

fn native_interval(
    prepared: &PreparedIncandescentSource<'_>,
    environment: &impl LumpedEndpointEnvironment,
    mut input: UniformEnthalpyStepInput<'_>,
    duration: f64,
    prescribed_power: f64,
    checkpoint: &mut impl FnMut() -> Result<(), ConductionError>,
) -> Result<IncandescentStep, IncandescentHistoryError> {
    let heat = prescribed_power * duration;
    if !duration.is_finite()
        || duration <= 0.0
        || !heat.is_finite()
        || (prescribed_power != 0.0 && heat == 0.0)
    {
        return Err(IncandescentHistoryError::Invalid(
            "native interval or prescribed heat is not representable",
        ));
    }
    input.duration_s = duration;
    input.internal_heat_j = heat;
    Ok(prepared.advance(environment, input, checkpoint)?)
}
