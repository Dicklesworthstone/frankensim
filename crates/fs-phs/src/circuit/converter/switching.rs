//! Explicit complementary switching on the original converter storage.
//!
//! Binary bridge connections use exactly the parent's original J/R/G and L/C
//! law. Each physical gate interval advances the complete native trajectory;
//! flux and bus charge are continuous at every commanded edge. The original
//! positive heat integral therefore includes switching ripple's conduction
//! loss and its correlation with each particular channel's conducting time.
//! No averaged-current square or fitted switching-loss map replaces that heat.
//!
//! This is a prescribed, ideal complementary gate schedule with finite linear
//! on-resistance, not a semiconductor transition model. It adds no dead time,
//! diode/recovery behavior, gate/output capacitance, switching-transition heat,
//! DCM clamp, automatic controller or device safe-area certificate. External
//! commands specify every edge; state-triggered event localization is absent.

use super::{
    AveragedBuckBoost, BridgeModulation, BuckBoostParameters, CONDUCTION_CHANNELS, ConductionHeat,
    ConverterStep, PreparedAveragedBuckBoost, invalid, poll, positive,
};
use crate::linear_exact::LinearSplitBudget;
use crate::{PhsError, PreparedStepError, PreparedStepRecord, storage_energy_change};

/// Bound the physical event work and native caches before preparation.
pub const MAX_SWITCH_INTERVALS: usize = 64;

/// One conducting device on each synchronous half bridge. Its complementary
/// device is OFF. Both devices off/on cannot be inferred from this state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BridgeSwitchState {
    pub input_high: bool,
    pub output_high: bool,
}

impl BridgeSwitchState {
    fn modulation(self) -> BridgeModulation {
        BridgeModulation {
            input_high_fraction: if self.input_high { 1.0 } else { 0.0 },
            output_high_fraction: if self.output_high { 1.0 } else { 0.0 },
        }
    }
}

/// A positive physical duration at one commanded connection. Intervals remain
/// in caller order, including repetitions that can carry different sources.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SwitchInterval {
    pub state: BridgeSwitchState,
    pub duration_s: f64,
}

/// Cold work/storage admission and the complete-cycle physical energy gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SwitchCycleControls {
    pub max_intervals: usize,
    pub absolute_energy_tolerance_j: f64,
    /// Relative to complete-cycle storage and absolute physical work/heat.
    /// This allowance is never increased to make a requested cycle pass.
    pub relative_energy_tolerance: f64,
}

impl Default for SwitchCycleControls {
    fn default() -> Self {
        Self {
            max_intervals: MAX_SWITCH_INTERVALS,
            absolute_energy_tolerance_j: 0.0,
            relative_energy_tolerance: 1e-12,
        }
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn total(self) -> f64 {
        self.value + self.correction
    }
}

/// Centered PWM on one common carrier: each positive-bus connection is high
/// on `[(1-d)*T/2, (1+d)*T/2]`. Simultaneous edges are grouped, with no invented
/// intermediate topology or zero-duration pulse. Exactly zero/one duty becomes
/// a static connection. The returned intervals preserve the represented clock;
/// [`PreparedSwitchCycle::period_s`] reports their compensated duration sum.
///
/// # Errors
/// Nonfinite/nonpositive period, duty outside [0,1], or a requested nonzero
/// high/low pulse whose edges cannot be represented distinctly in f64.
pub fn centered_pwm(
    period_s: f64,
    modulation: BridgeModulation,
) -> Result<Vec<SwitchInterval>, PhsError> {
    let duty = [
        modulation.input_high_fraction,
        modulation.output_high_fraction,
    ];
    if !positive(period_s)
        || duty
            .iter()
            .any(|d| !d.is_finite() || !(0.0..=1.0).contains(d))
    {
        return Err(invalid(
            "centered PWM requires positive period and duties in [0,1]",
        ));
    }
    let mut windows = [[0.0; 2]; 2];
    let mut edges = vec![0.0, period_s];
    for side in 0..2 {
        if duty[side] == 1.0 {
            windows[side] = [0.0, period_s];
        } else if duty[side] > 0.0 {
            let start = (0.5 * (1.0 - duty[side])) * period_s;
            let end = period_s - start;
            if !(start > 0.0 && start < end && end < period_s) {
                return Err(invalid("centered PWM pulse edges are not representable"));
            }
            windows[side] = [start, end];
            edges.extend([start, end]);
        }
    }
    edges.sort_by(f64::total_cmp);
    edges.dedup();
    let mut intervals = Vec::with_capacity(edges.len() - 1);
    for edge in edges.windows(2) {
        let high = |side: usize| {
            duty[side] > 0.0 && edge[0] >= windows[side][0] && edge[1] <= windows[side][1]
        };
        let duration_s = edge[1] - edge[0];
        if !positive(duration_s) {
            return Err(invalid("centered PWM interval is unrepresentable"));
        }
        intervals.push(SwitchInterval {
            state: BridgeSwitchState {
                input_high: high(0),
                output_high: high(1),
            },
            duration_s,
        });
    }
    Ok(intervals)
}

struct CachedInterval {
    interval: SwitchInterval,
    native: PreparedAveragedBuckBoost,
}

/// A reusable physical gate cycle. Native kernels with identical connection and
/// duration share one cache. Every authored interval and input slot survives;
/// no original electrical state is projected or reset at a gate transition.
///
/// Preparation allocates and builds exponentials. Successful warm calls perform
/// O(intervals) native three-state updates without new allocation/exponentials.
/// This is an operation-count statement, not a measured real-time guarantee.
pub struct PreparedSwitchCycle {
    parameters: BuckBoostParameters,
    intervals: Vec<SwitchInterval>,
    cache_indices: Vec<usize>,
    caches: Vec<CachedInterval>,
    period_s: f64,
    averaged_modulation: BridgeModulation,
    controls: SwitchCycleControls,
}

impl PreparedSwitchCycle {
    /// Prepare every distinct original fixed-connection interval. The unchanged
    /// native scalar-state budget applies to each kernel; `max_intervals` bounds
    /// the number of kernels and all per-cycle event work before allocation.
    /// Material L/C and all resistances are held over the whole cycle. Changing
    /// them between cycles requires the owning physical material transaction.
    pub fn new_controlled(
        parameters: BuckBoostParameters,
        intervals: &[SwitchInterval],
        budget: LinearSplitBudget,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        if !(1..=MAX_SWITCH_INTERVALS).contains(&controls.max_intervals)
            || intervals.is_empty()
            || intervals.len() > controls.max_intervals
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !controls.relative_energy_tolerance.is_finite()
            || controls.relative_energy_tolerance < 0.0
        {
            return Err(invalid("invalid switching cycle interval or energy allowance").into());
        }
        let mut duration = Sum::default();
        let mut high_time = [Sum::default(); 2];
        for interval in intervals {
            poll(&mut cancelled)?;
            if !positive(interval.duration_s) {
                return Err(invalid("switch intervals require positive finite duration").into());
            }
            duration.add(interval.duration_s);
            if interval.state.input_high {
                high_time[0].add(interval.duration_s);
            }
            if interval.state.output_high {
                high_time[1].add(interval.duration_s);
            }
        }
        let period_s = duration.total();
        if !positive(period_s) {
            return Err(invalid("switching cycle period is unrepresentable").into());
        }
        let averaged_modulation = BridgeModulation {
            input_high_fraction: high_time[0].total() / period_s,
            output_high_fraction: high_time[1].total() / period_s,
        };
        let mut caches: Vec<CachedInterval> = Vec::new();
        let mut cache_indices = Vec::with_capacity(intervals.len());
        for &interval in intervals {
            poll(&mut cancelled)?;
            let cached = caches.iter().position(|cache| {
                cache.interval.state == interval.state
                    && cache.interval.duration_s.to_bits() == interval.duration_s.to_bits()
            });
            let index = match cached {
                Some(index) => index,
                None => {
                    let native = AveragedBuckBoost::new(parameters, interval.state.modulation())?
                        .prepare_controlled(interval.duration_s, budget, &mut cancelled)?;
                    let index = caches.len();
                    caches.push(CachedInterval { interval, native });
                    index
                }
            };
            cache_indices.push(index);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            parameters,
            intervals: intervals.to_vec(),
            cache_indices,
            caches,
            period_s,
            averaged_modulation,
            controls,
        })
    }

    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    pub fn intervals(&self) -> &[SwitchInterval] {
        &self.intervals
    }
    pub const fn period_s(&self) -> f64 {
        self.period_s
    }
    /// Physical connection-time fractions, useful for an independently admitted
    /// averaged comparison. This does not claim the two trajectories coincide.
    pub const fn averaged_modulation(&self) -> BridgeModulation {
        self.averaged_modulation
    }
    pub fn cached_interval_count(&self) -> usize {
        self.caches.len()
    }

    /// One full cycle with the same held terminal currents throughout. Returned
    /// voltages are actual cycle means and work is the full trajectory integral.
    /// Neither output is published if any interval or final gate refuses.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        terminal_current_a: &[f64; 2],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConverterStep, PreparedStepError> {
        self.advance(
            initial,
            |_| *terminal_current_a,
            next,
            average_bus_voltage_v,
            cancelled,
        )
    }

    /// The caller may instead hold a DIFFERENT current pair on each original
    /// gate interval. Each source multiplies its own interval voltage integral;
    /// averaging sources before computing work would discard their correlation.
    /// The exact interval count/order is available through [`Self::intervals`].
    pub fn step_interval_sources_controlled(
        &mut self,
        initial: &[f64; 3],
        terminal_current_a: &[[f64; 2]],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConverterStep, PreparedStepError> {
        if terminal_current_a.len() != self.intervals.len() {
            return Err(
                invalid("one current pair is required per original switch interval").into(),
            );
        }
        self.advance(
            initial,
            |interval| terminal_current_a[interval],
            next,
            average_bus_voltage_v,
            cancelled,
        )
    }

    fn advance(
        &mut self,
        initial: &[f64; 3],
        current: impl Fn(usize) -> [f64; 2],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConverterStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let (mut heat, mut voltage_integral, mut work) = (
            [Sum::default(); CONDUCTION_CHANNELS],
            [Sum::default(); 2],
            [Sum::default(); 2],
        );
        let (mut native_heat, mut native_work, mut native_change, mut absolute_work) = (
            Sum::default(),
            Sum::default(),
            Sum::default(),
            Sum::default(),
        );
        let mut storage_scale = self.caches[0].native.system().hamiltonian(initial).abs();
        let mut last_voltage = [0.0; 2];
        for k in 0..self.intervals.len() {
            poll(&mut cancelled)?;
            let (mut endpoint, mut voltage) = ([0.0; 3], [0.0; 2]);
            let record = self.caches[self.cache_indices[k]].native.step_controlled(
                &state,
                &current(k),
                &mut endpoint,
                &mut voltage,
                &mut cancelled,
            )?;
            for channel in 0..CONDUCTION_CHANNELS {
                heat[channel].add(record.conduction.joules[channel]);
            }
            for port in 0..2 {
                voltage_integral[port].add(self.intervals[k].duration_s * voltage[port]);
                work[port].add(record.terminal_work_j[port]);
                absolute_work.add(record.terminal_work_j[port].abs());
            }
            native_heat.add(record.native.dissipated);
            native_work.add(record.native.supplied);
            native_change.add(record.native.delta_h);
            state = endpoint;
            last_voltage = voltage;
            storage_scale =
                storage_scale.max(self.caches[0].native.system().hamiltonian(&state).abs());
        }
        let system = self.caches[0].native.system();
        let delta_h = storage_energy_change(
            system.storage(),
            initial,
            &state,
            &mut [0.0; 3],
            &mut [0.0; 3],
        );
        let joules = heat.map(Sum::total);
        let terminal_work_j = work.map(Sum::total);
        // A static one-interval cycle retains the native mean exactly instead
        // of multiplying and dividing it by the same rounded duration again.
        let voltage = if self.intervals.len() == 1 {
            last_voltage
        } else {
            voltage_integral.map(|value| value.total() / self.period_s)
        };
        let dissipated = native_heat.total();
        let supplied = native_work.total();
        let scale = storage_scale
            .max(absolute_work.total())
            .max(dissipated.abs())
            .max(delta_h.abs());
        let tolerance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        let residuals = [
            delta_h + dissipated - supplied,
            delta_h - native_change.total(),
            joules.iter().sum::<f64>() - dissipated,
            terminal_work_j.iter().sum::<f64>() - supplied,
        ];
        if !scale.is_finite()
            || !tolerance.is_finite()
            || voltage
                .iter()
                .chain(&joules)
                .chain(&terminal_work_j)
                .any(|value| !value.is_finite())
            || joules.iter().any(|value| *value < 0.0)
            || residuals
                .iter()
                .any(|value| !value.is_finite() || value.abs() > tolerance)
        {
            return Err(
                invalid("complete switching cycle physical work/heat balance refused").into(),
            );
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&state);
        average_bus_voltage_v.copy_from_slice(&voltage);
        Ok(ConverterStep {
            native: PreparedStepRecord {
                delta_h,
                dissipated,
                supplied,
                newton_iters: 0,
                solver_residual: 0.0,
                solver_certificate: None,
            },
            terminal_work_j,
            conduction: ConductionHeat { joules },
        })
    }
}
