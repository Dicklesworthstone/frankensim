//! Prescribed switching with simultaneous source impedance and load response.
//!
//! Each interval uses the original voltage-driven converter, including both
//! bus capacitors and the inductor. Source/load currents respond throughout the
//! interval; sampling them once at an edge would change the physical circuit.
//! Winding/channel heat and both terminal-resistor heats integrate their own
//! squared currents on that same trajectory. An ideal source can receive
//! regenerated energy, while a zero-volt source is an ordinary resistive load.
//!
//! Gates remain ideal complementary commands. This composes the existing
//! finite linear conduction model, not a semiconductor transition, diode,
//! dead-time, electrochemical storage or temperature-evolution model.

use super::{
    AveragedBuckBoost, BridgeModulation, BuckBoostParameters, CONDUCTION_CHANNELS, ConductionHeat,
    ConverterStep, Sum, SwitchCycleControls, SwitchInterval, admit_cycle, invalid, poll,
};
use crate::circuit::converter::terminals::{PreparedVoltageDrivenBuckBoost, VoltageDrivenStep};
use crate::linear_exact::{DenseExponentialBudget, LinearSplitBudget};
use crate::{PreparedStepError, PreparedStepRecord, storage_energy_change};

struct CachedInterval {
    interval: SwitchInterval,
    native: PreparedVoltageDrivenBuckBoost,
}

/// Cached switching cycle with a finite Thevenin source/load on each bus.
///
/// All material parameters and terminal resistances are held over the cycle.
/// Original flux and charges survive every edge without projection. Identical
/// connection/duration pairs share a propagator and heat observers. Successful
/// warm steps allocate nothing and require O(intervals) fixed-size updates;
/// source values and interval order remain distinct even when caches are shared.
pub struct PreparedVoltageDrivenSwitchCycle {
    parameters: BuckBoostParameters,
    terminal_resistance_ohm: [f64; 2],
    intervals: Vec<SwitchInterval>,
    cache_indices: Vec<usize>,
    caches: Vec<CachedInterval>,
    period_s: f64,
    averaged_modulation: BridgeModulation,
    controls: SwitchCycleControls,
}

impl PreparedVoltageDrivenSwitchCycle {
    /// Prepare the existing simultaneous loaded trajectory for each distinct
    /// gate interval. The native and heat budgets apply per distinct cache;
    /// the cycle interval cap bounds their count before allocation.
    ///
    /// # Errors
    /// Invalid cycle/physical parameters, unavailable native work or storage
    /// budget, unrepresentable arithmetic, or cancellation. No partial owner
    /// is returned. Semiconductor operating limits remain caller supplied.
    pub fn new_controlled(
        parameters: BuckBoostParameters,
        terminal_resistance_ohm: [f64; 2],
        intervals: &[SwitchInterval],
        budget: LinearSplitBudget,
        heat_budget: DenseExponentialBudget,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        let (period_s, averaged_modulation) = admit_cycle(intervals, controls, &mut cancelled)?;
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
                        .prepare_voltage_driven_controlled(
                        terminal_resistance_ohm,
                        interval.duration_s,
                        budget,
                        heat_budget,
                        &mut cancelled,
                    )?;
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
            terminal_resistance_ohm,
            intervals: intervals.to_vec(),
            cache_indices,
            caches,
            period_s,
            averaged_modulation,
            controls,
        })
    }

    /// Frozen original storage and conduction parameters.
    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    /// Physical source/load resistance on each bus [ohm].
    pub const fn terminal_resistance_ohm(&self) -> [f64; 2] {
        self.terminal_resistance_ohm
    }
    /// All original gate intervals, in authored order.
    pub fn intervals(&self) -> &[SwitchInterval] {
        &self.intervals
    }
    /// Compensated sum of the represented gate durations [s].
    pub const fn period_s(&self) -> f64 {
        self.period_s
    }
    /// Connection-time fractions for an independently admitted averaged model.
    /// Equality of fractions does not establish equal ripple or equal losses.
    pub const fn averaged_modulation(&self) -> BridgeModulation {
        self.averaged_modulation
    }
    /// Number of distinct native connection/duration caches.
    pub fn cached_interval_count(&self) -> usize {
        self.caches.len()
    }

    /// Advance a cycle under two held ideal source voltages. Terminal currents
    /// evolve with both buses throughout every interval; they are never frozen
    /// to an edge sample. Output voltages are complete-cycle time means.
    /// Cancellation/refusal preserves both caller buffers and permits retry.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        source_voltage_v: &[f64; 2],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<VoltageDrivenStep, PreparedStepError> {
        self.advance(
            initial,
            |_| *source_voltage_v,
            next,
            average_bus_voltage_v,
            cancelled,
        )
    }

    /// Hold a separately supplied voltage pair on each original gate interval.
    /// Source work integrates each voltage with its own current response,
    /// retaining gate/source correlation and signed regenerative work.
    pub fn step_interval_sources_controlled(
        &mut self,
        initial: &[f64; 3],
        source_voltage_v: &[[f64; 2]],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<VoltageDrivenStep, PreparedStepError> {
        if source_voltage_v.len() != self.intervals.len() {
            return Err(
                invalid("one voltage pair is required per original switch interval").into(),
            );
        }
        self.advance(
            initial,
            |interval| source_voltage_v[interval],
            next,
            average_bus_voltage_v,
            cancelled,
        )
    }

    fn advance(
        &mut self,
        initial: &[f64; 3],
        source: impl Fn(usize) -> [f64; 2],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VoltageDrivenStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let mut heat = [Sum::default(); CONDUCTION_CHANNELS];
        let mut terminal_heat = [Sum::default(); 2];
        let mut source_work = [Sum::default(); 2];
        let mut terminal_work = [Sum::default(); 2];
        let mut voltage_integral = [Sum::default(); 2];
        let mut native_change = Sum::default();
        let mut native_heat = Sum::default();
        let mut native_work = Sum::default();
        let mut absolute_work = Sum::default();
        let mut absolute_change = Sum::default();
        let mut last_voltage = [0.0; 2];
        for k in 0..self.intervals.len() {
            poll(&mut cancelled)?;
            let (mut endpoint, mut voltage) = ([0.0; 3], [0.0; 2]);
            let record = self.caches[self.cache_indices[k]].native.step_controlled(
                &state,
                &source(k),
                &mut endpoint,
                &mut voltage,
                &mut cancelled,
            )?;
            for channel in 0..CONDUCTION_CHANNELS {
                heat[channel].add(record.converter.conduction.joules[channel]);
            }
            for port in 0..2 {
                terminal_heat[port].add(record.terminal_resistor_heat_j[port]);
                source_work[port].add(record.source_work_j[port]);
                terminal_work[port].add(record.converter.terminal_work_j[port]);
                voltage_integral[port].add(self.intervals[k].duration_s * voltage[port]);
                absolute_work.add(record.source_work_j[port].abs());
                absolute_work.add(record.converter.terminal_work_j[port].abs());
            }
            native_change.add(record.converter.native.delta_h);
            absolute_change.add(record.converter.native.delta_h.abs());
            native_heat.add(record.converter.native.dissipated);
            native_work.add(record.converter.native.supplied);
            state = endpoint;
            last_voltage = voltage;
        }
        let delta_h = storage_energy_change(
            self.caches[0].native.system().storage(),
            initial,
            &state,
            &mut [0.0; 3],
            &mut [0.0; 3],
        );
        let joules = heat.map(Sum::total);
        let terminal_resistor_heat_j = terminal_heat.map(Sum::total);
        let source_work_j = source_work.map(Sum::total);
        let terminal_work_j = terminal_work.map(Sum::total);
        let voltage = if self.intervals.len() == 1 {
            last_voltage
        } else {
            voltage_integral.map(|value| value.total() / self.period_s)
        };
        let dissipated = native_heat.total();
        let supplied = native_work.total();
        let source_total = source_work_j.iter().sum::<f64>();
        let terminal_heat_total = terminal_resistor_heat_j.iter().sum::<f64>();
        // Use physical transfers in THIS cycle, never a large unchanged
        // capacitor energy datum, to set the relative error allowance.
        let scale = absolute_work
            .total()
            .max(absolute_change.total())
            .max(dissipated.abs())
            .max(terminal_heat_total)
            .max(delta_h.abs());
        let tolerance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        let residuals = [
            delta_h + dissipated - supplied,
            delta_h + dissipated + terminal_heat_total - source_total,
            delta_h - native_change.total(),
            joules.iter().sum::<f64>() - dissipated,
            terminal_work_j.iter().sum::<f64>() - supplied,
        ];
        if !scale.is_finite()
            || !tolerance.is_finite()
            || voltage
                .iter()
                .chain(&joules)
                .chain(&terminal_resistor_heat_j)
                .chain(&source_work_j)
                .chain(&terminal_work_j)
                .any(|v| !v.is_finite())
            || joules
                .iter()
                .chain(&terminal_resistor_heat_j)
                .any(|v| *v < 0.0)
            || residuals
                .iter()
                .any(|v| !v.is_finite() || v.abs() > tolerance)
        {
            return Err(invalid(
                "complete voltage-driven switching cycle physical work/heat balance refused",
            )
            .into());
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&state);
        average_bus_voltage_v.copy_from_slice(&voltage);
        Ok(VoltageDrivenStep {
            converter: ConverterStep {
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
            },
            source_work_j,
            terminal_resistor_heat_j,
        })
    }
}
