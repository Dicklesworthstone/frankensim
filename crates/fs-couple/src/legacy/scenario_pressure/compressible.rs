//! Time-dependent pressure supplies on the original conservative gas clock.
//!
//! A stationary external reservoir supplies an explicit pressure and temperature
//! history in the field's original frozen composition. The native characteristic
//! flux selects inflow, outflow and choking. Its total enthalpy, including
//! pressure work, is the sole boundary energy transfer; no separate p dV term,
//! nominal-air density or incompressible pressure-to-density mapping is added.
//!
//! Every numerical retry samples its own SSP-RK2 stage times. Steps stop at the
//! next authored signal breakpoint. A held-table jump is approached from the
//! left at the final stage and from the right on the following initial stage,
//! so a held pressure pulse cannot be skipped or injected half a step early.
//! Smooth histories still require temporal refinement; these are prescribed
//! infinite supplies, not finite gas storage or a nonreflecting termination.

use fs_gas::dynamics::{
    Conserved, Error as GasError, Euler1d, FaceFlux, FrozenGas, PortStage, StepControl, StepReport,
};
use fs_qty::{Dims, Pressure, Temperature, Time};
use fs_scenario::signal::DEFAULT_MAXIMUM_CYCLES_PER_STEP;
use fs_scenario::{Interp, ScenarioError, TimeSignal};
use std::cell::RefCell;

/// Explicit interpretation of the pressure history before source evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ReservoirPressureReference {
    /// Values are positive absolute pressure [Pa].
    Absolute,
    /// Values are signed gauge pressure [Pa] above this positive absolute
    /// reference. Their sum, not the gauge value, must remain positive.
    Gauge { reference_pressure_pa: f64 },
}

/// Borrowed histories of one stationary, infinite external gas reservoir.
///
/// Pressure carries Pa dimensions and temperature carries K dimensions.
/// During inflow these are stagnation conditions; during outflow the source
/// supplies backpressure while the field supplies the outgoing entropy and
/// enthalpy. The composition and energy reference are explicitly shared with
/// the field. No exhaled mixture, moving reservoir or implicit temperature is
/// inferred. Array index zero binds the left end and index one the right end.
#[derive(Clone, Copy, Debug)]
pub struct ReservoirHistory<'a> {
    pub pressure: &'a TimeSignal,
    pub pressure_reference: ReservoirPressureReference,
    pub temperature: &'a TimeSignal,
}

/// Actual source and original native face flux used in one accepted RK stage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirStage {
    /// Physical stage time [s]. At a held-table jump the final stage uses the
    /// interval's left limit; the next initial stage uses the new held value.
    pub time_s: f64,
    pub absolute_pressure_pa: f64,
    pub prescribed_temperature_k: f64,
    /// Native per-unit-area flux along the increasing duct coordinate.
    /// At the left end positive mass enters; at the right end positive leaves.
    pub flux: FaceFlux,
}

/// Only accepted stage samples are returned; rejected trial samples disappear.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirHistoryStep {
    /// Native integrated exchanges and conservation residual, including all
    /// original unselected boundaries. Energy is in the source's native gauge.
    pub gas: StepReport,
    /// `[RK stage][left or right end]`; `None` retains the original boundary.
    pub stages: [[Option<ReservoirStage>; 2]; 2],
}

/// A refusal preserves the original field, time and boundary configuration.
#[derive(Debug)]
pub enum ReservoirHistoryError {
    Invalid(&'static str),
    Signal(ScenarioError),
    Gas(GasError),
    #[cfg(feature = "nonlinear-brass")]
    Duct(fs_duct::nonlinear::NonlinearDuctError),
    Cancelled,
}

impl core::fmt::Display for ReservoirHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "compressible pressure history: {self:?}")
    }
}
impl std::error::Error for ReservoirHistoryError {}
impl From<GasError> for ReservoirHistoryError {
    fn from(error: GasError) -> Self {
        match error {
            GasError::Cancelled => Self::Cancelled,
            error => Self::Gas(error),
        }
    }
}

/// Advance an original Euler interval with prescribed pressure/temperature.
///
/// The borrowed source must equal the field's complete source, including its
/// caloric domain and formation-energy reference. Retain that immutable source
/// when constructing the field; this adapter never clones gas cards or cell
/// storage in the step loop. The original `TimeSignal` evaluator owns units,
/// validation, interpolation and clamped tails. Input-history/source refusals
/// terminate the call; numerical face/CFL refusals retain native retries.
///
/// Controls and unsupported periodic/hydrostatic end ports retain the native
/// gas admission. The stored boundary configuration is not rebound on success:
/// subsequent driven steps must continue through this function.
pub fn step_euler(
    field: &mut Euler1d,
    source: &FrozenGas,
    histories: [Option<ReservoirHistory<'_>>; 2],
    maximum_duration: Time,
    controls: StepControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ReservoirHistoryStep, ReservoirHistoryError> {
    poll(&mut cancelled)?;
    if source != field.gas() {
        return Err(ReservoirHistoryError::Invalid(
            "reservoir histories require the field's complete original gas source",
        ));
    }
    advance(
        source,
        histories,
        field.time().value(),
        maximum_duration,
        |duration, ports, check| {
            field
                .step_with_ports(duration, controls, ports, check)
                .map_err(ReservoirHistoryError::from)
        },
        cancelled,
    )
}

/// The same stage-synchronous exchange through the retained nonlinear duct
/// geometry and cell-observation API. Enabled with `nonlinear-brass`, the
/// existing feature that admits this crate's nonlinear Euler duct dependency.
#[cfg(feature = "nonlinear-brass")]
pub fn step_duct(
    duct: &mut fs_duct::nonlinear::NonlinearDuct,
    source: &FrozenGas,
    histories: [Option<ReservoirHistory<'_>>; 2],
    maximum_duration: Time,
    controls: StepControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ReservoirHistoryStep, ReservoirHistoryError> {
    use fs_duct::nonlinear::NonlinearDuctError;
    poll(&mut cancelled)?;
    if source != duct.solver().gas() {
        return Err(ReservoirHistoryError::Invalid(
            "reservoir histories require the duct's complete original gas source",
        ));
    }
    advance(
        source,
        histories,
        duct.time().value(),
        maximum_duration,
        |duration, ports, check| {
            duct.step_with_ports(duration, controls, ports, check)
                .map_err(|error| match error {
                    NonlinearDuctError::Gas(error) => error.into(),
                    error => ReservoirHistoryError::Duct(error),
                })
        },
        cancelled,
    )
}

type Ports<'a> = dyn FnMut(&PortStage) -> Result<[Option<FaceFlux>; 2], GasError> + 'a;

fn advance(
    source: &FrozenGas,
    histories: [Option<ReservoirHistory<'_>>; 2],
    time_s: f64,
    maximum_duration: Time,
    native: impl FnOnce(
        Time,
        &mut Ports<'_>,
        &mut dyn FnMut() -> bool,
    ) -> Result<StepReport, ReservoirHistoryError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ReservoirHistoryStep, ReservoirHistoryError> {
    if histories.iter().all(Option::is_none)
        || !maximum_duration.value().is_finite()
        || maximum_duration.value() <= 0.0
    {
        return Err(ReservoirHistoryError::Invalid(
            "pressure history needs a selected reservoir and positive finite duration",
        ));
    }
    let mut duration = maximum_duration.value();
    for history in histories.iter().flatten() {
        for signal in [history.pressure, history.temperature] {
            if let Some(knot) = next_breakpoint(signal, time_s) {
                duration = duration.min(knot - time_s);
            }
            if let TimeSignal::Chirp(chirp) = signal {
                // The native solver may shorten this interval further for CFL
                // or retry, but its RK stage samples must first resolve every
                // selected pressure and temperature source's physical frequency.
                let end = chirp.resolved_step_end(
                    time_s,
                    time_s + duration,
                    DEFAULT_MAXIMUM_CYCLES_PER_STEP,
                ).map_err(ReservoirHistoryError::Signal)?;
                duration = duration.min(end - time_s);
            }
        }
    }
    let cancelled = RefCell::new(cancelled);
    let failure = RefCell::new(None);
    let mut stages = [[None; 2]; 2];
    let result = {
        let mut ports = |stage: &PortStage| {
            if stage.stage == 0 {
                stages = [[None; 2]; 2];
            }
            let mut fluxes = [None; 2];
            let mut accepted = [None; 2];
            for (end, history) in histories.iter().enumerate() {
                let Some(history) = history else { continue };
                let mut check = || (cancelled.borrow_mut())();
                let (reservoir, pressure, temperature) =
                    match sample(source, *history, stage.time, stage.stage == 1, &mut check) {
                        Ok(value) => value,
                        Err(error) => {
                            // A malformed prescribed history cannot be repaired by
                            // halving the numerical timestep. Retain its diagnostic
                            // while the native cancellation path aborts atomically.
                            *failure.borrow_mut() = Some(error);
                            return Err(GasError::Cancelled);
                        }
                    };
                let (interior, normal) = if end == 0 {
                    (stage.left.0, -1.0)
                } else {
                    (stage.right.0, 1.0)
                };
                let flux = source.reservoir_flux(interior, reservoir, normal, &mut check)?;
                fluxes[end] = Some(flux);
                accepted[end] = Some(ReservoirStage {
                    time_s: stage.time,
                    absolute_pressure_pa: pressure,
                    prescribed_temperature_k: temperature,
                    flux,
                });
            }
            stages[stage.stage] = accepted;
            Ok(fluxes)
        };
        let mut check = || failure.borrow().is_some() || (cancelled.borrow_mut())();
        native(Time::new(duration), &mut ports, &mut check)
    };
    if let Some(error) = failure.into_inner() {
        return Err(error);
    }
    Ok(ReservoirHistoryStep {
        gas: result?,
        stages,
    })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReservoirHistoryError> {
    if cancelled() {
        Err(ReservoirHistoryError::Cancelled)
    } else {
        Ok(())
    }
}

fn sample(
    source: &FrozenGas,
    history: ReservoirHistory<'_>,
    time_s: f64,
    left_limit: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Conserved, f64, f64), ReservoirHistoryError> {
    let supplied = value(
        history.pressure,
        Pressure::DIMS,
        time_s,
        left_limit,
        cancelled,
    )?;
    let pressure = match history.pressure_reference {
        ReservoirPressureReference::Absolute => supplied,
        ReservoirPressureReference::Gauge {
            reference_pressure_pa,
        } => {
            if !reference_pressure_pa.is_finite() || reference_pressure_pa <= 0.0 {
                return Err(ReservoirHistoryError::Invalid(
                    "gauge reservoir pressure needs a positive absolute reference",
                ));
            }
            reference_pressure_pa + supplied
        }
    };
    let temperature = value(
        history.temperature,
        Temperature::DIMS,
        time_s,
        left_limit,
        cancelled,
    )?;
    if !pressure.is_finite() || pressure <= 0.0 || !temperature.is_finite() || temperature <= 0.0 {
        return Err(ReservoirHistoryError::Invalid(
            "reservoir history requires positive absolute pressure and temperature",
        ));
    }
    poll(cancelled)?;
    let density = (pressure / source.gas_constant()) / temperature;
    let reservoir = source.conserved(density, 0.0, pressure)?;
    poll(cancelled)?;
    Ok((reservoir, pressure, temperature))
}

fn value(
    signal: &TimeSignal,
    dimensions: Dims,
    time_s: f64,
    left_limit: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ReservoirHistoryError> {
    let mut evaluation_time = time_s;
    if left_limit
        && let TimeSignal::Table {
            times,
            interp: Interp::Hold,
            ..
        } = signal
    {
        let index = times.partition_point(|time| *time < time_s);
        if index > 0 && times.get(index) == Some(&time_s) {
            evaluation_time = time_s.next_down();
        }
    }
    let quantity = signal
        .eval_with_checkpoint(evaluation_time, &mut |_| poll(cancelled))?
        .map_err(ReservoirHistoryError::Signal)?;
    if quantity.dims != dimensions {
        return Err(ReservoirHistoryError::Invalid(
            "reservoir history must carry pressure in Pa and temperature in K",
        ));
    }
    Ok(quantity.value)
}

// Scheduling reads the same immutable signal; its complete semantic validation
// still occurs in the original evaluator before any gas state can be accepted.
fn next_breakpoint(signal: &TimeSignal, time_s: f64) -> Option<f64> {
    let next = match signal {
        TimeSignal::Constant(_) => None,
        TimeSignal::Chirp(chirp) => [chirp.t_start_s(), chirp.t_end_s()]
            .into_iter()
            .find(|&knot| knot > time_s),
        TimeSignal::Ramp { t_start, t_end, .. } => {
            [*t_start, *t_end].into_iter().find(|time| *time > time_s)
        }
        TimeSignal::Table { times, .. } => times
            .get(times.partition_point(|time| *time <= time_s))
            .copied(),
        TimeSignal::Chebfun(profile) => {
            let (start, end) = profile.cheb.domain();
            [start, end].into_iter().find(|time| *time > time_s)
        }
    };
    next.filter(|time| time.is_finite() && *time > time_s)
}

#[cfg(test)]
mod tests;
