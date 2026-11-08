//! Held electrical drive histories on the original implicit thermal owner.
//!
//! `TimeSignal` owns the SI schedule and right-continuous hold convention.
//! Each incoming interval retains its held drive through its endpoint; a jump
//! changes the next interval's drive without resetting the body's enthalpy.
//! The explicit maximum step controls temporal resolution independently of
//! native source-feedback admission and energy-residual tolerances.

use super::{IncandescentError, IncandescentStep, PreparedIncandescentSource};
use fs_conduction::ConductionError;
use fs_conduction::lumped::{LumpedThermalEnvironment, MAX_LUMPED_ENTHALPY_STEPS};
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

/// Physical clock and bounded thermal work; no error controller is implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncandescentHistoryControls {
    /// Absolute start in the original signal's seconds-based clock.
    pub start_s: f64,
    /// Positive maximum physical substep [s], also clipped at every drive knot.
    pub maximum_step_s: f64,
    /// Bounds native solves and retained rows. A held table may contain at
    /// most this many entries plus one, including entries outside the horizon.
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
    environment: &LumpedThermalEnvironment,
    input: UniformEnthalpyStepInput<'_>,
    drive_history: IncandescentDriveHistory<'_>,
    controls: IncandescentHistoryControls,
    mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
) -> Result<Vec<TimedIncandescentStep>, IncandescentHistoryError> {
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
    let mut samples = Vec::new();
    while time < end {
        checkpoint()?;
        if samples.len() == controls.maximum_steps {
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
        while time < span_end {
            checkpoint()?;
            if samples.len() == controls.maximum_steps {
                return Err(Invalid(
                    "incandescent history exhausted its native step budget",
                ));
            }
            let remaining = span_end - time;
            let maximum_duration = limit.min(remaining);
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
            let heat = prescribed_power * duration;
            if !next.is_finite()
                || !duration.is_finite()
                || duration <= 0.0
                || !heat.is_finite()
                || (prescribed_power != 0.0 && heat == 0.0)
            {
                return Err(Invalid(
                    "native substep or prescribed heat cannot be represented on this clock",
                ));
            }
            samples
                .try_reserve(1)
                .map_err(|_| Invalid("history output allocation failed"))?;
            current.duration_s = duration;
            current.internal_heat_j = heat;
            let step = prepared.advance(environment, current, &mut checkpoint)?;
            current.initial = step.thermal.state;
            samples.push(TimedIncandescentStep {
                start_s: time,
                end_s: next,
                step,
            });
            time = next;
        }
        if time < end {
            held = values[next_event];
            next_event += 1;
        }
    }
    checkpoint()?;
    Ok(samples)
}
